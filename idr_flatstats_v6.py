from omero.gateway import BlitzGateway
from PIL import Image as PILImage
import numpy as np
import matplotlib.pyplot as plt
import math

# === CONNECTION ===
conn = BlitzGateway('public', 'public',
                    host='idr.openmicroscopy.org',
                    port=4064, secure=True)
conn.connect()

# === 1. PARAMETERS ===
image_id = 2858397   # change to the desired IDR image ID

img    = conn.getObject("Image", image_id)
pixels = img.getPrimaryPixels()

pixel_size_xy = img.getPixelSizeX()   # µm
pixel_size_z  = img.getPixelSizeZ()   # µm

if pixel_size_xy is None:
    raise ValueError("pixel_size_xy missing")
if pixel_size_z is None or pixel_size_z == 0:
    raise ValueError("pixel_size_z missing or zero")

size_x = img.getSizeX()
size_y = img.getSizeY()
size_z = img.getSizeZ()
print(f"Image: {img.getName()} | {size_x}x{size_y}, Z={size_z}, "
      f"pxy={pixel_size_xy} µm, pz={pixel_size_z} µm")

# === 2. CENTRAL Z-SLICE SELECTION ===
# Load a central subset of planes covering at least 2.5 µm
Z_THICKNESS_UM = 2.5
Z_MIN_PLANES   = 5
Z_MAX_PLANES   = 50

n_planes  = min(Z_MAX_PLANES, max(Z_MIN_PLANES, math.ceil(Z_THICKNESS_UM / pixel_size_z)))
z_center  = size_z // 2
z_start   = max(0, z_center - n_planes // 2)
z_end     = min(size_z, z_start + n_planes)
z_indices = list(range(z_start, z_end))
print(f"Z slice: {z_start}–{z_end-1} ({len(z_indices)} planes out of {size_z})")

# === 3. STACK LOADING (channel 0) ===
stack = np.array([pixels.getPlane(z, 0, 0) for z in z_indices])
print(f"Stack loaded: {stack.shape}, dtype={stack.dtype}")

# === 4. SCALING ===
# Downsample to ~100 px along the short axis
target = 100
smaller_dim = min(size_x, size_y)
scaling_factor = target / smaller_dim
if scaling_factor > 1.0:
    raise ValueError(f"Image too small ({smaller_dim}px < {target}px)")

print(f"Scaling factor: {scaling_factor:.4f}")

nz, ny, nx = stack.shape
new_ny = int(ny * scaling_factor)
new_nx = int(nx * scaling_factor)
stack_scaled = np.array([
    np.array(PILImage.fromarray(stack[z].astype(np.float32)).resize(
        (new_nx, new_ny), PILImage.Resampling.BILINEAR))
    for z in range(nz)
])
pixel_size_xy_scaled = pixel_size_xy / scaling_factor

# === 5. Z-MAP (for display) ===
# For each pixel, record the z position of maximum intensity
max_proj = stack_scaled.max(axis=0)
best_z   = np.argmax(stack_scaled, axis=0)
z_map    = best_z * pixel_size_z

# Masked z-map for display (top 10% of max projection)
threshold_display = np.percentile(max_proj, 90)
mask_display      = max_proj > threshold_display
z_map_masked      = np.where(mask_display, z_map, 0)

# === 6. ROBUSTNESS ANALYSIS (percentile thresholds 88 to 92) ===
# Plane fitting repeated at five percentile thresholds for robustness
results = []   # list of (slope, dir_trigo) per threshold

for p in range(88, 93):
    threshold = np.percentile(max_proj, p)
    mask_p    = max_proj > threshold
    if mask_p.sum() < 10:
        continue
    y_c, x_c = np.where(mask_p)
    x_p = x_c * pixel_size_xy_scaled
    y_p = y_c * pixel_size_xy_scaled
    z_p = z_map[y_c, x_c]
    mx, my, mz = x_p.mean(), y_p.mean(), z_p.mean()
    dx, dy, dz = x_p-mx, y_p-my, z_p-mz
    sxx=(dx*dx).sum(); syy=(dy*dy).sum(); sxy=(dx*dy).sum()
    sxz=(dx*dz).sum(); syz=(dy*dz).sum()
    det = sxx*syy - sxy*sxy
    if det == 0:
        continue
    a = (syy*sxz - sxy*syz) / det
    b = (sxx*syz - sxy*sxz) / det
    slope     = np.sqrt(a**2+b**2) * 100
    dir_trigo = np.degrees(np.arctan2(b, a))
    results.append((slope, dir_trigo))

# === 7. SUMMARY ===
slope_mean = np.mean([r[0] for r in results])
slope_sd   = np.std([r[0] for r in results])
x_vals  = [slope * np.cos(np.radians(dir_t)) for slope, dir_t in results]
y_vals  = [slope * np.sin(np.radians(dir_t)) for slope, dir_t in results]
dir_ref = np.degrees(np.arctan2(np.mean(y_vals), np.mean(x_vals)))
diffs   = [((a - dir_ref + 180) % 360) - 180 for a in [r[1] for r in results]]
dir_sd  = np.sqrt(np.mean([d**2 for d in diffs]))
dir_compass = (90 + dir_ref) % 360   # convert trigonometric to compass bearing

print(f"=== RESULTS ===")
print(f"Slope     : {slope_mean:.2f} ± {slope_sd:.2f} µm/100µm")
print(f"Direction : {dir_compass:.1f} ± {dir_sd:.1f}°")

# === 8. OUTPUT FILENAME ===
project    = img.getProject()
project_id = project.getId() if project else "noproject"
image_title = img.getName().replace(" ", "_").replace("/", "_")
filename = f"{project_id}_{image_id}_{image_title}.png"

# === 9. VISUALIZATION ===
ny_s, nx_s = z_map_masked.shape
fig, ax = plt.subplots(figsize=(6, 6))
im = ax.imshow(z_map_masked, cmap='inferno', origin='upper')
plt.colorbar(im, ax=ax, label='Z position (µm)')

# Individual percentile estimates (white crosses)
for slope, dir_t in results:
    angle_rad = np.radians(dir_t)
    ax.plot(nx_s/2 + 25*np.cos(angle_rad),
            ny_s/2 + 25*np.sin(angle_rad),
            'w+', markersize=12, markeredgewidth=2)

# Mean tilt vector (white arrow)
angle_rad = np.radians(dir_ref)
length = min(nx_s, ny_s) / 4
cx, cy = nx_s/2, ny_s/2
ax.annotate("",
    xy=(cx + np.cos(angle_rad)*length, cy + np.sin(angle_rad)*length),
    xytext=(cx - np.cos(angle_rad)*length, cy - np.sin(angle_rad)*length),
    arrowprops=dict(arrowstyle='->', color='white', lw=3, mutation_scale=40))

ax.set_title(f"{img.getName()}\n"
             f"Slope: {slope_mean:.2f}±{slope_sd:.2f} µm/100µm | "
             f"Dir: {dir_compass:.1f}±{dir_sd:.1f}°")
plt.tight_layout()
plt.savefig(filename, dpi=150, bbox_inches='tight')
plt.show()
print(f"Saved: {filename}")
conn.close()
