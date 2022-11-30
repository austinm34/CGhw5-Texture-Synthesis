# CGhw5
extra hw

All source is in CGhw5/CGhw5

I used lodepng.h to load pngs into a vector

From there I stored them int 2d vectors for easier operation

I use a 3x3 seed to start synthesis

I iterate throught every pixel where the window will still be in bounds of the sample image when trying to find matching neighborhoods to the pixel to be synthesized.

This causes the algorithm to be extremly, unfathomably slow.

Right now only window size is 3x3. it was made with generalization in mind and this can easily be added to where window size is free.

Could not get program to finish in time to upload an output.

The only output image is sand_out.png which is the seed as well a single pixel that was synthesized as a test. it is below the seed in the image.
