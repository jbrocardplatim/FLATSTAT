from omero.gateway import BlitzGateway
import csv
import os

# =============================================================================
# PARAMETERS
# =============================================================================
OUTPUT_CSV = "idr_census.csv"          # Output census file (adapt path if needed)
LOG_FILE   = "idr_census_done.txt"     # Ratchet file (processed project IDs)
MAX_IMAGES_PER_PROJECT = 10000         # Skip projects exceeding this image count

# Imaging method filters
IMAGING_METHODS_OK      = ['confocal', 'fluorescence']
IMAGING_METHODS_EXCLUDE = ['light sheet', 'sheet', 'lightsheet', 'spim', 'lsfm']

# Free-text annotation keys to search for instrument information
FREETEXT_KEYS = ['protocol', 'comment', 'experiment comment',
                 'method', 'microscope', 'acquisition']

FIELDNAMES = [
    'project_id', 'project_name', 'release_date', 'imaging_method', 'pubmed_id',
    'n_datasets', 'n_images_total', 'n_images_eligible',
    'dataset_id', 'dataset_name',
    'image_id', 'image_name',
    'size_x', 'size_y', 'size_z', 'size_t', 'size_c',
    'pixel_size_xy', 'pixel_size_z',
    'obj_magnification', 'obj_na', 'obj_model', 'obj_manufacturer',
    'freetext_instrument'
]

# =============================================================================
# RATCHET: SKIP ALREADY PROCESSED PROJECTS
# =============================================================================
done_project_ids = set()
if os.path.exists(LOG_FILE):
    with open(LOG_FILE, 'r') as f:
        for line in f:
            line = line.strip()
            if line:
                done_project_ids.add(int(line))
print(f"Resuming — {len(done_project_ids)} projects already processed")

# =============================================================================
# OMERO CONNECTION
# =============================================================================
conn = BlitzGateway('public', 'public',
                    host='idr.openmicroscopy.org',
                    port=4064, secure=True)

try:
    conn.connect()
    print("Connected to IDR")

    projects = list(conn.getObjects("Project"))
    print(f"{len(projects)} projects found")

    file_exists = os.path.exists(OUTPUT_CSV)
    csvfile = open(OUTPUT_CSV, 'a', newline='', encoding='utf-8')
    writer  = csv.DictWriter(csvfile, fieldnames=FIELDNAMES)
    if not file_exists:
        writer.writeheader()

    for project in projects:

        proj_id = project.getId()

        # Skip already processed projects
        if proj_id in done_project_ids:
            print(f"[SKIP ratchet] {project.getName()}")
            continue

        # === Project-level metadata ===
        imaging_method = ''
        release_date   = ''
        pubmed_id      = ''

        for ann in project.listAnnotations():
            if hasattr(ann, 'getValue'):
                pairs = ann.getValue()
                if isinstance(pairs, list):
                    for key, val in pairs:
                        kl = key.lower()
                        if kl == 'imaging method':
                            imaging_method = val
                        elif kl == 'release date':
                            release_date = val
                        elif kl == 'pubmed id':
                            pubmed_id = val

        # Imaging method filter
        im_lower = imaging_method.lower()
        if any(ex in im_lower for ex in IMAGING_METHODS_EXCLUDE):
            print(f"[SKIP excluded] {project.getName()} — {imaging_method}")
            with open(LOG_FILE, 'a') as f:
                f.write(f"{proj_id}\n")
            continue

        if not any(m in im_lower for m in IMAGING_METHODS_OK):
            print(f"[SKIP method]  {project.getName()} — {imaging_method or 'absent'}")
            with open(LOG_FILE, 'a') as f:
                f.write(f"{proj_id}\n")
            continue

        print(f"\n[OK] {project.getName()} | {imaging_method} | {release_date}")

        # === PASS 1: count images and detect oversized projects ===
        datasets       = list(project.listChildren())
        n_datasets     = len(datasets)
        n_images_total = 0
        oversized      = False

        for dataset in datasets:
            images = list(dataset.listChildren())
            n_images_total += len(images)
            if n_images_total > MAX_IMAGES_PER_PROJECT:
                oversized = True
                break

        if oversized:
            print(f"  -> [SKIP oversized] {n_images_total}+ images, project skipped")
            with open(LOG_FILE, 'a') as f:
                f.write(f"{proj_id}\n")
            continue

        # === PASS 2: write CSV rows (up to 5 eligible images per dataset) ===
        n_eligible = 0
        for dataset in datasets:
            n_images_dataset = 0
            for img in dataset.listChildren():
                if n_images_dataset >= 5:
                    break

                sx = img.getSizeX()
                sy = img.getSizeY()
                sz = img.getSizeZ()

                # Eligibility filters
                if sz <= 5:
                    continue
                if min(sx, sy) <= 100:
                    continue
                pxy = img.getPixelSizeX()
                if pxy is None:
                    continue

                n_eligible       += 1
                n_images_dataset += 1

                # Instrument metadata — source 1: structured OMERO instrument object
                obj_mag, obj_na, obj_model, obj_manuf = '', '', '', ''
                instr = img.getInstrument()
                if instr:
                    for obj in instr.getObjectives():
                        obj_mag   = obj.getNominalMagnification() or ''
                        obj_na    = obj.getLensNA() or ''
                        obj_model = obj.getModel() or ''
                        obj_manuf = obj.getManufacturer() or ''
                        break

                # Instrument metadata — source 2: free-text annotations
                freetext_instrument = ''
                for ann in img.listAnnotations():
                    if hasattr(ann, 'getValue'):
                        pairs = ann.getValue()
                        if isinstance(pairs, list):
                            for key, val in pairs:
                                if any(k in key.lower() for k in FREETEXT_KEYS):
                                    freetext_instrument += f"{key}: {val} | "

                row = {
                    'project_id':          proj_id,
                    'project_name':        project.getName(),
                    'release_date':        release_date,
                    'imaging_method':      imaging_method,
                    'pubmed_id':           pubmed_id,
                    'n_datasets':          n_datasets,
                    'n_images_total':      n_images_total,
                    'n_images_eligible':   '',
                    'dataset_id':          dataset.getId(),
                    'dataset_name':        dataset.getName(),
                    'image_id':            img.getId(),
                    'image_name':          img.getName(),
                    'size_x':              sx,
                    'size_y':              sy,
                    'size_z':              sz,
                    'size_t':              img.getSizeT(),
                    'size_c':              img.getSizeC(),
                    'pixel_size_xy':       pxy,
                    'pixel_size_z':        img.getPixelSizeZ() or '',
                    'obj_magnification':   obj_mag,
                    'obj_na':              obj_na,
                    'obj_model':           obj_model,
                    'obj_manufacturer':    obj_manuf,
                    'freetext_instrument': freetext_instrument.strip(' |')
                }
                writer.writerow(row)
                csvfile.flush()

        # Mark project as processed
        with open(LOG_FILE, 'a') as f:
            f.write(f"{proj_id}\n")

        print(f"  -> {n_datasets} datasets | {n_images_total} images | "
              f"{n_eligible} eligible")

finally:
    csvfile.close()
    conn.close()
    print("\nDone.")
