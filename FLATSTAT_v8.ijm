macro "FLATSTAT Action Tool - C000T4b12F" {

	// === CHECKS ===
	if (nImages == 0 || nSlices < 3) exit("No Z-stack open. Please get a Z-stack (nSlices >= 3) ready."); 
	Stack.getDimensions(width, height, channels, slices, frames);
	if (channels > 1) run("Split Channels");

	// === INITIALIZATION ===
	dir = getDirectory("image");
	print(dir);
	if (dir=="") dir = getDirectory("Choose output directory");
	getVoxelSize(pixelWidth, pixelHeight, zstep, unit);
	print("=== FLATSTAT ANALYSIS ===");

	// === 1. PRE-PROCESSING ===
	saveAs("Tiff",dir+getTitle());
	t = getTitle();

	// Compute scaling factor to reach ~100 pixels along the short axis
	w0 = getWidth();
	h0 = getHeight();
	short_dim = minOf(w0, h0);
	if (short_dim < 100) exit("Image too small: shortest dimension (" + short_dim + " px) must be >= 100 px.");
	scaling_factor = 100 / short_dim;

	print("Original size: " + w0 + " x " + h0 + " px");
	print("Scaling factor: " + d2s(scaling_factor, 4));

	run("Scale...", "x=" + scaling_factor + " y=" + scaling_factor + " z=1.0 interpolation=Bilinear average process create");
	rename("scaled");

	// Retrieve scaled dimensions
	w = getWidth();
	h = getHeight();
	nZ = nSlices;
	pixelWidth_scaled = pixelWidth / scaling_factor;

	// === 2. Z-MAP COMPUTATION WITH MASK ===
	// Build mask from top ~10% brightest pixels of max-intensity projection
	selectWindow("scaled");
	run("Z Project...", "projection=[Max Intensity]");
	rename("MAX_scaled");
	run("Enhance Contrast...", "saturated=10");
	getMinAndMax(min, max);
	if (max==255) setMinAndMax(min, 1);
	run("Apply LUT");
	run("8-bit");

	// Create Z-map
	newImage("Z_Map", "32-bit", w, h, 1);

	count_valid = 0;
	max_points  = w * h;
	x_list = newArray(max_points);
	y_list = newArray(max_points);
	z_list = newArray(max_points);

	print("Computing Z-map...");
	setBatchMode(true);
	for (y = 0; y < h; y++) {
		if (y % 10 == 0) print("  row " + y + "/" + h);
		for (x = 0; x < w; x++) {
			selectWindow("MAX_scaled");
			mask_value = getPixel(x, y);
			if (mask_value > 200) {
				best_intensity = 0;
				best_z = 0;
				selectWindow("scaled");
				for (z = 1; z <= nZ; z++) {
					setSlice(z);
					intensity = getPixel(x, y);
					if (intensity > best_intensity) {
						best_intensity = intensity;
						best_z = z;
					}
				}
				selectWindow("Z_Map");
				z_position = best_z * zstep;
				setPixel(x, y, z_position);
				x_list[count_valid] = x * pixelWidth_scaled;
				y_list[count_valid] = y * pixelWidth_scaled;
				z_list[count_valid] = z_position;
				count_valid++;
			}
		}
	}
	setBatchMode(false);
	print("Valid pixels: " + count_valid);

	// Display Z-map
	selectWindow("Z_Map");
	run("Fire");
	resetMinAndMax();

	// === 3. PLANE FITTING AND SLOPE COMPUTATION ===
	if (count_valid > 100) {
		print("Fitting best-fit plane...");

		// Centroid
		mx = 0; my = 0; mz = 0;
		for (i = 0; i < count_valid; i++) {
			mx += x_list[i]; my += y_list[i]; mz += z_list[i];
		}
		mx /= count_valid; my /= count_valid; mz /= count_valid;

		// Least-squares coefficients: z = a*x + b*y + c (determinant method)
		sxx = 0; syy = 0; sxy = 0; sxz = 0; syz = 0;
		for (i = 0; i < count_valid; i++) {
			dx = x_list[i] - mx; dy = y_list[i] - my; dz = z_list[i] - mz;
			sxx += dx*dx; syy += dy*dy; sxy += dx*dy;
			sxz += dx*dz; syz += dy*dz;
		}

		det = sxx*syy - sxy*sxy;
		if (det != 0) {
			a = (syy*sxz - sxy*syz) / det;
			b = (sxx*syz - sxy*sxz) / det;

			// Slope magnitude and compass direction
			pente          = sqrt(a*a + b*b);
			slope_100um = pente * 100;
			direction      = atan2(b, a) * 180 / PI;
			dir_compass_ij = (90 + direction + 360) % 360;

			print("=== RESULTS ===");
			print("Slope: " + d2s(slope_100um, 2) + " µm/100µm");
			print("Direction (compass): " + d2s(dir_compass_ij, 1) + "°");

			// Rescale Z-map to ~1000 px along the largest dimension
			selectWindow("Z_Map");
			target_size   = 1000;
			max_dimension = maxOf(w, h);
			scale_back    = target_size / max_dimension;
			run("Scale...", "x=" + scale_back + " y=" + scale_back + " interpolation=Bilinear create title=" + substring(t, 0, lengthOf(t)-4) + "_zmap.tif");

			// Draw tilt direction arrow
			img_w = getWidth(); img_h = getHeight();
			cx = img_w / 2;    cy = img_h / 2;
			angleRad = dir_compass_ij * PI / 180;
			dx = sin(angleRad); dy = -cos(angleRad);
			L  = 100 * floor(cx / 100);
			makeArrow(cx, cy, cx + L*dx, cy + L*dy, "filled");
			Roi.setStrokeWidth(4);
			run("Add Selection...");

		} else {
			exit("ERROR: Singular matrix — plane fit impossible.");
		}
	} else {
		exit("ERROR: Not enough valid pixels...");
	}

	// === 4. ANNOTATION AND SAVE ===
	setFont("Arial", 22, "Bold");
	setColor(0, 255, 255);
	drawString("Slope: "     + d2s(slope_100um, 2) + " µm/100µm", 24, 24);
	drawString("Direction: " + d2s(dir_compass_ij, 1) + "°",         24, 48);

	tz = getTitle();
	saveAs("Tiff", dir + tz);

	selectWindow("Log");
	saveAs("Text", dir + substring(tz, 0, lengthOf(tz)-4) + ".txt");
	run("Close");

	// Cleanup
	close("Z_Map"); close("MAX_scaled"); close("scaled");
}
