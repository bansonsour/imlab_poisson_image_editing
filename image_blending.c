#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <dirent.h>
#include "imcore.h"
#include "cvcore.h"

// struct to keep four neighbours of the (i,j) and itself
struct sparse_element {
    int x;
    int y;
    double value[5];
};

void sparse_element_set(struct sparse_element *out, int x, int y, int idx, double value)
{
    out->x = x;
    out->y = y;
    out->value[idx] = value;
}

matrix_t* poisson_solver(vector_t *possion_matrix, matrix_t *gradient)
{
    int dx[5] = {0, -1, +1, 0, 0};
    int dy[5] = {0, 0, 0, -1, +1};

    double w = 1.8;
    int i, j, k, t, Max_Iter = 4000;

    // get the address of the data
    struct sparse_element *poisson_values = vdata(possion_matrix, 0);

    // allocate memory for the solution
    matrix_t *solution = matrix_create(gradient, NULL);
    
    double gradient_value;

    for (k = 0; k < Max_Iter; k++)
    {
        for (i = 0; i < length(possion_matrix); i++)
        {
            double s = 0;
            for (j = 1; j < 5; j++)
            {
                if (poisson_values[i].value[j] > 0)
                {
                    s += poisson_values[i].value[j] * at(double, solution, poisson_values[i].y + dy[j], poisson_values[i].x + dx[j]);
                }  
            }

            // get the gradient value @(y,x)
            matrix_get(gradient, poisson_values[i].y, poisson_values[i].x, 0, &gradient_value);

            double val = at(double, solution, poisson_values[i].y, poisson_values[i].x);
            val = (1 - w) * val + (w / poisson_values[i].value[0]) * (gradient_value + s);
            matrix_set(solution, poisson_values[i].y, poisson_values[i].x, 0, &val);
        }
    }

    return solution;
}

// blend the given two images using the poisson image editing algortihm
return_t poisson_blend(matrix_t *background, matrix_t *foreground, matrix_t *mask, int x1, int y1, int mode, matrix_t *destination)
{
    int i = 0, j = 0, c = 0, k = 0;

    // copy the data of the background into the destination
    matrix_copy(background, destination);

    int dx[5] = {0, -1, +1, 0, 0};
    int dy[5] = {0, 0, 0, -1, +1};
   
    // solve Poisson equation for each channel
    for (c = 0; c < channels(background); c++)
    {
        // create 1d gradient vectors
        matrix_t *gradient = matrix_create(double, rows(foreground), cols(foreground), 1);

        // allocate possion matrix (sparse format)
        vector_t *possion_matrix = vector_create(struct sparse_element);

        for (i = 1; i < height(foreground) - 1; i++)
        {
            for (j = 1; j < width(foreground) - 1; j++)
            {
                // if (i,j) is inside the mask
                if(at(uint8_t, mask, i,j,0) > 0)
                {
                    struct sparse_element poisson_value = {0};
                    double gradient_value = 0;

                    sparse_element_set(&poisson_value, j, i, 0, 4);

                    // compute the average gradient at(i,j) using the four corner of the pixel
                    for (k = 1; k < 5; k++)
                    {
                        // if we are on border, force I(x,y) = T(x,y)
                        //if (y == 0 || y == height(foreground) - 1 || x == 0 || x == width(foreground) - 1)
                        if (!(at(uint8_t, mask, i + dy[k], j + dx[k], 0) > 0))
                        {
                            gradient_value += (double)at(uint8_t, background, i + dy[k] + y1, j + dx[k] + x1, c);
                        }
                        else
                        {
                            sparse_element_set(&poisson_value, j, i, k, 1);
                        }

                        // mix the gradient of the background and foreground
                        double f1 = (double)at(uint8_t, foreground, i, j, c) - (double)at(uint8_t, foreground, i + dy[k], j + dx[k], c);
                        double g1 = (double)at(uint8_t, background, i + y1, j + x1, c) - (double)at(uint8_t, background, i + dy[k] + y1, j + dx[k] + x1, c);

                        if(mode == 0)
                        {
                            gradient_value += f1;
                        }
                        else
                        {
                            gradient_value += fabs(f1) > fabs(g1) ? f1 : g1;
                        }
                        
                    }

                    matrix_set(gradient, i, j, 0, &gradient_value);
                    vector_push(possion_matrix, &poisson_value);
                }
            }
        }
        // solve poisson equation
        matrix_t *solution = poisson_solver(possion_matrix, gradient);

        struct sparse_element *poisson_values = vdata(possion_matrix, 0);

        int t = 0;
        for (t = 0; t < length(possion_matrix); t++)
        {
            double val;
            matrix_get(solution, poisson_values[t].y, poisson_values[t].x, 0, &val);
            at(uint8_t, destination, poisson_values[t].y + y1, poisson_values[t].x + x1, c) = clamp(val, 0, 255);
        }

        // free vector and possion matrix
        matrix_free(&gradient);
        matrix_free(&solution);
        vector_free(&possion_matrix);
    }

    // return success
    return SUCCESS;
}

int main(int argc, char *argv[]) 
{
    if(argc != 5)
    {
        printf("call image_blending folder x y mode\n");
        return -1;
    }

    // get the file names
    char background_filename[512];
    char foreground_filename[512];
    char foreground_mask_filename[512];

    snprintf(background_filename, 512, "%s//background.bmp", argv[1]);
    snprintf(foreground_filename, 512, "%s//foreground.bmp", argv[1]);
    snprintf(foreground_mask_filename, 512, "%s//foreground_mask.bmp", argv[1]);

    // x and y position of the foreground
    int x = atoi(argv[2]);
    int y = atoi(argv[3]);
    int mode = atoi(argv[4]);

    // read foreground and background images
    matrix_t *background = imread(background_filename);
    matrix_t *foreground = imread(foreground_filename);
    matrix_t *foreground_mask = imread(foreground_mask_filename);

    if (background == NULL || foreground == NULL || foreground_mask  == NULL)
    {
        printf("unable to open input files\n");
        return -1;
    }

    // create an empty matrix
    matrix_t *blended = matrix_create(uint8_t);

    // run the blending algorithm
    poisson_blend(background, foreground, foreground_mask, x, y, mode, blended);

    // write the output image
    imwrite(blended, "blended_image.bmp");

    matrix_free(&background);
    matrix_free(&foreground);
    matrix_free(&foreground_mask);
    matrix_free(&blended);

    return 0;
}