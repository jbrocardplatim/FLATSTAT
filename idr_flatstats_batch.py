from omero.gateway import BlitzGateway
from PIL import Image as PILImage
import numpy as np
import csv
import math
import os

# =============================================================================
# PARAMETERS
# =============================================================================
INPUT_CSV   = "idr_census.csv"           # Census CSV (adapt path if needed)
OUTPUT_CSV  = "flatstat_results.csv"     # FlatStat results output
LOG_FILE    = "flatstat_done.txt"        # Ratchet file 

SAVE_FIGURES = False                     # True = save one figure per stack
if SAVE_FIGURES:
    import matplotlib.pyplot as plt
FIGURES_DIR  = "flatstat_figures"        # Output folder for figures

# FlatStat eligibility filters (redundant with census, applied here as safeguard
# and for positive imaging_method filter)
IMAGING_METHOD_OK = ['confocal', 'fluorescence']
MAX_PIXEL_SIZE_XY = 10.0   # µm
MAX_PIXEL_SIZE_Z  = 10.0   # µm
MAX_SIZE_XY       = 4096   # px (max side)

# Central z-slice parameters
Z_THICKNESS_UM = 2.5       # target thickness in µm
Z_MIN_PLANES   = 5         # minimum number of planes
Z_MAX_PLANES   = 50        # maximum number of planes

FIELDNAMES = [
    'image_id', 'image_name', 'project_name', 'dataset_id',
    'channel', 'n_planes_loaded',
    'pixel_size_xy', 'pixel_size_z',
    'slope_mean', 'slope_sd',
    'dir_compass', 'dir_sd',
    'status'
]

# =============================================================================
# FUNCTIONS
# =============================================================================

def compute_n_planes(pixel_size_z):
    """Number of planes to load based on z-step size."""
    n = math.ceil(Z_THICKNESS_UM / pixel_size_z)
    return min(Z_MAX_PLANES, max(Z_MIN_PLANES, n))


def flatstat(stack, pixel_size_xy, pixel_size_z):
    """
    Compute slope and tilt direction from a numpy stack (nz, ny, nx).
    Returns (slope_mean, slope_sd, dir_compass, dir_sd, z_map, max_proj, results, dir_ref)
    or raises ValueError.
    """
    target = 100
    nz, ny, nx = stack.shape
    smaller_dim = min(nx, ny)
    scaling_factor = target / smaller_dim

    # Bilinear downscaling
    new_ny = int(ny * scaling_factor)
    new_nx = int(nx * scaling_factor)
    stack_scaled = np.array([
        np.array(PILImage.fromarray(stack[z].astype(np.float32)).resize(
            (new_nx, new_ny), PILImage.Resampling.BILINEAR))
        for z in range(nz)
    ])
    pixel_size_xy_scaled = pixel_size_xy / scaling_factor

    # Z-map: for each pixel, z position of maximum intensity
    max_proj = stack_scaled.max(axis=0)
    best_z   = np.argmax(stack_scaled, axis=0)
    z_map    = best_z * pixel_size_z

    # Robustness: repeat plane fitting at percentile thresholds 88-92
    results = []
    for p in range(88, 93):
        threshold = np.percentile(max_proj, p)
        mask_p = max_proj > threshold
        if mask_p.sum() < 10:
            continue
        y_c, x_c = np.where(mask_p)
        x_p = x_c * pixel_size_xy_scaled
        y_p = y_c * pixel_size_xy_scaled
        z_p = z_map[y_c, x_c]
        mx, my, mz = x_p.mean(), y_p.mean(), z_p.mean()
        dx, dy, dz = x_p - mx, y_p - my, z_p - mz
        sxx = (dx*dx).sum(); syy = (dy*dy).sum(); sxy = (dx*dy).sum()
        sxz = (dx*dz).sum(); syz = (dy*dz).sum()
        det = sxx*syy - sxy*sxy
        if det == 0:
            continue
        a = (syy*sxz - sxy*syz) / det
        b = (sxx*syz - sxy*sxz) / det
        slope     = np.sqrt(a**2 + b**2) * 100
        dir_trigo = np.degrees(np.arctan2(b, a))
        results.append((slope, dir_trigo))

    if len(results) < 3:
        raise ValueError(f"Only {len(results)} valid percentile thresholds")

    # Aggregate results across percentile thresholds
    slope_mean = np.mean([r[0] for r in results])
    slope_sd   = np.std([r[0] for r in results])
    x_vals = [r[0] * np.cos(np.radians(r[1])) for r in results]
    y_vals = [r[0] * np.sin(np.radians(r[1])) for r in results]
    dir_ref     = np.degrees(np.arctan2(np.mean(y_vals), np.mean(x_vals)))
    diffs       = [((a - dir_ref + 180) % 360) - 180 for a in [r[1] for r in results]]
    dir_sd      = np.sqrt(np.mean([d**2 for d in diffs]))
    dir_compass = (90 + dir_ref) % 360   # convert trigonometric to compass bearing

    return slope_mean, slope_sd, dir_compass, dir_sd, z_map, max_proj, results, dir_ref


def save_figure(image_name, z_map_masked, results, dir_ref,
                slope_mean, slope_sd, dir_compass, dir_sd, channel, out_dir):
    """Save a diagnostic figure showing the Z-map and tilt vector."""
    ny_s, nx_s = z_map_masked.shape
    fig, ax = plt.subplots(figsize=(6, 6))
    im = ax.imshow(z_map_masked, cmap='inferno', origin='upper')
    plt.colorbar(im, ax=ax, label='Z position (µm)')
    # Individual percentile estimates (crosses)
    for slope, dir_t in results:
        angle_rad = np.radians(dir_t)
        ax.plot(nx_s/2 + 25*np.cos(angle_rad),
                ny_s/2 + 25*np.sin(angle_rad),
                'w+', markersize=12, markeredgewidth=2)
    # Mean tilt vector (arrow)
    angle_rad = np.radians(dir_ref)
    length = min(nx_s, ny_s) / 4
    cx, cy = nx_s/2, ny_s/2
    ax.annotate("",
        xy=(cx + np.cos(angle_rad)*length, cy + np.sin(angle_rad)*length),
        xytext=(cx - np.cos(angle_rad)*length, cy - np.sin(angle_rad)*length),
        arrowprops=dict(arrowstyle='->', color='white', lw=3, mutation_scale=40))
    ax.set_title(f"{image_name} | c={channel}\n"
                 f"Slope: {slope_mean:.2f}±{slope_sd:.2f} µm/100µm | "
                 f"Dir: {dir_compass:.1f}±{dir_sd:.1f}°")
    plt.tight_layout()
    safe_name = image_name.replace('/', '_').replace(' ', '_')
    path = os.path.join(out_dir, f"{safe_name}_c{channel}.png")
    plt.savefig(path, dpi=80)
    plt.close()


# =============================================================================
# LOAD AND FILTER CENSUS CSV
# =============================================================================
def safe_float(val):
    try:
        return float(val)
    except (ValueError, TypeError):
        return None

eligible = []
seen_ids = set()

with open(INPUT_CSV, newline='', encoding='utf-8') as f:
    reader = csv.DictReader(f)
    for row in reader:
        # Normalize freetext column name (may be corrupted by Excel)
        for key in list(row.keys()):
            if key.startswith('freetext') and key != 'freetext_instrument':
                row['freetext_instrument'] = row.pop(key)

        # Positive imaging method filter
        method = (row.get('imaging_method') or '').lower()
        if not any(m in method for m in IMAGING_METHOD_OK):
            continue

        # Numeric filters
        pxy = safe_float(row.get('pixel_size_xy'))
        pz  = safe_float(row.get('pixel_size_z'))
        sx  = safe_float(row.get('size_x'))
        sy  = safe_float(row.get('size_y'))
        iid = safe_float(row.get('image_id'))

        if pxy is None or pxy > MAX_PIXEL_SIZE_XY:
            continue
        if pz is None or pz <= 0 or pz > MAX_PIXEL_SIZE_Z:
            continue
        if sx is None or sy is None or max(sx, sy) > MAX_SIZE_XY:
            continue
        if iid is None:
            continue

        row['image_id']      = int(iid)
        row['dataset_id']    = int(safe_float(row.get('dataset_id')) or 0)
        row['pixel_size_xy'] = pxy
        row['pixel_size_z']  = pz
        eligible.append(row)
        seen_ids.add(int(iid))

print(f"Eligible images: {len(eligible)} ({len(seen_ids)} unique image_ids)")

# =============================================================================
# RATCHET: SKIP ALREADY PROCESSED IMAGES
# =============================================================================
done_ids = set()
if os.path.exists(LOG_FILE):
    with open(LOG_FILE) as f:
        for line in f:
            line = line.strip()
            if line:
                done_ids.add(int(line))
print(f"Resuming — {len(done_ids)} images already processed")

# =============================================================================
# OMERO CONNECTION
# =============================================================================
conn = BlitzGateway('public', 'public',
                    host='idr.openmicroscopy.org',
                    port=4064, secure=True)
conn.connect()
print("Connected to IDR")

if SAVE_FIGURES:
    os.makedirs(FIGURES_DIR, exist_ok=True)

# =============================================================================
# MAIN BATCH LOOP
# =============================================================================
file_exists = os.path.exists(OUTPUT_CSV)
csvfile = open(OUTPUT_CSV, 'a', newline='', encoding='utf-8')
writer  = csv.DictWriter(csvfile, fieldnames=FIELDNAMES)
if not file_exists:
    writer.writeheader()

try:
    for row in eligible:
        image_id     = int(row['image_id'])
        project_name = row['project_name']
        dataset_id   = int(row['dataset_id'])

        # Skip already processed images
        if image_id in done_ids:
            print(f"  [SKIP] {image_id}")
            continue

        print(f"\n[->] image_id={image_id} | {project_name}")

        try:
            img    = conn.getObject("Image", image_id)
            pixels = img.getPrimaryPixels()

            pxy  = img.getPixelSizeX()
            pz   = img.getPixelSizeZ()
            sx   = img.getSizeX()
            sy   = img.getSizeY()
            sz   = img.getSizeZ()
            sc   = img.getSizeC()
            name = img.getName()

            # Consistency checks (values may differ from CSV)
            if pxy is None or pz is None or pz == 0:
                raise ValueError("Missing pixel size from OMERO metadata")
            if max(sx, sy) > MAX_SIZE_XY:
                raise ValueError(f"Image too large: {sx}x{sy}")

            # Central z-slice
            n_planes  = compute_n_planes(pz)
            z_center  = sz // 2
            z_start   = max(0, z_center - n_planes // 2)
            z_end     = min(sz, z_start + n_planes)
            z_indices = list(range(z_start, z_end))
            print(f"    {name} | {sx}x{sy} Z={sz} C={sc} "
                  f"| pxy={pxy:.4f} pz={pz:.4f} | z slice={z_start}-{z_end-1} ({len(z_indices)} planes)")

            # Loop over channels
            for c in range(sc):
                try:
                    stack = np.array([pixels.getPlane(z, c, 0) for z in z_indices])

                    slope_mean, slope_sd, dir_compass, dir_sd, \
                        z_map, max_proj, results, dir_ref = flatstat(stack, pxy, pz)

                    print(f"    c={c} -> slope={slope_mean:.2f}±{slope_sd:.2f} "
                          f"µm/100µm | dir={dir_compass:.1f}±{dir_sd:.1f}°")

                    if SAVE_FIGURES:
                        threshold_display = np.percentile(max_proj, 90)
                        z_map_masked = np.where(max_proj > threshold_display, z_map, 0)
                        save_figure(name, z_map_masked, results, dir_ref,
                                    slope_mean, slope_sd, dir_compass, dir_sd,
                                    c, FIGURES_DIR)

                    writer.writerow({
                        'image_id':        image_id,
                        'image_name':      name,
                        'project_name':    project_name,
                        'dataset_id':      dataset_id,
                        'channel':         c,
                        'n_planes_loaded': len(z_indices),
                        'pixel_size_xy':   round(pxy, 6),
                        'pixel_size_z':    round(pz,  6),
                        'slope_mean':      round(slope_mean, 4),
                        'slope_sd':        round(slope_sd,   4),
                        'dir_compass':     round(dir_compass, 2),
                        'dir_sd':          round(dir_sd,      2),
                        'status':          'ok'
                    })

                except Exception as e_chan:
                    print(f"    c={c} -> ERROR: {e_chan}")
                    writer.writerow({
                        'image_id':        image_id,
                        'image_name':      name,
                        'project_name':    project_name,
                        'dataset_id':      dataset_id,
                        'channel':         c,
                        'n_planes_loaded': len(z_indices),
                        'pixel_size_xy':   round(pxy, 6),
                        'pixel_size_z':    round(pz,  6),
                        'slope_mean':      '',
                        'slope_sd':        '',
                        'dir_compass':     '',
                        'dir_sd':          '',
                        'status':          f'error: {e_chan}'
                    })

        except Exception as e_img:
            print(f"  ERROR image {image_id}: {e_img}")
            writer.writerow({
                'image_id':        image_id,
                'image_name':      row.get('image_name') or '',
                'project_name':    project_name,
                'dataset_id':      dataset_id,
                'channel':         '',
                'n_planes_loaded': '',
                'pixel_size_xy':   '',
                'pixel_size_z':    '',
                'slope_mean':      '',
                'slope_sd':        '',
                'dir_compass':     '',
                'dir_sd':          '',
                'status':          f'error: {e_img}'
            })

        # Mark image as processed (success or error)
        csvfile.flush()
        with open(LOG_FILE, 'a') as f:
            f.write(f"{image_id}\n")
        done_ids.add(image_id)

finally:
    csvfile.close()
    conn.close()
    print("\nDone.")
