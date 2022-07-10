#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

typedef struct{
	double x,y,z;
}vector;

void get_trajectory(double length, double width, double depth, int num_bodies, double body_radius, double body_mass, double time_step, vector pos[], vector vel[], vector force[], int n_sims){
    int i,j,k;
    double dist, mag, wtime;
    vector r, vhs, temp;
    FILE *ptr;
    ptr = fopen("trajectory.bin", "wb");
    if (ptr == NULL){
        printf("File not created!");
    }
    wtime = omp_get_wtime();
     
    for(i=0; i<n_sims; i++){
    #pragma omp parallel shared(n_sims, num_bodies, body_mass, body_radius, pos, vel, force, i) private(j,k,r,dist,mag, vhs) num_threads(8)
    {
        #pragma omp for collapse(2)
        for(j=0; j<num_bodies; j++){
            for(k=0; k<num_bodies; k++){
                if((i+1)%100 == 0 && j == 0 && k == 0)
                    fprintf(ptr, "%s %d\n", "Coordinates at time step", i+1);
                
                if (j != k){
                    r.x = pos[k].x - pos[j].x;
                    r.y = pos[k].y - pos[j].y;
                    r.z = pos[k].z - pos[j].z;
                    dist = sqrt(r.x*r.x + r.y*r.y + r.z*r.z);
                    if (dist == 0)
                        continue;
                    mag = (body_mass*body_mass)/(dist*dist);
                    force[j].x += (mag*r.x/dist);
                    force[j].y += (mag*r.y/dist);
                    force[j].z += (mag*r.z/dist);
                }   
            }

        }
        #pragma omp barrier

        #pragma omp for 
        for(j=0; j<num_bodies; j++){
            vhs.x = vel[j].x + (force[j].x*time_step)/(2*body_mass);
            vhs.y = vel[j].y + (force[j].y*time_step)/(2*body_mass);
            vhs.z = vel[j].z + (force[j].z*time_step)/(2*body_mass);

            pos[j].x += (vhs.x*time_step);
            pos[j].y += (vhs.y*time_step);
            pos[j].z += (vhs.z*time_step);

            if ((pos[j].x - body_radius) < 0)
                pos[j].x = body_radius;
            else if ((pos[j].x + body_radius) > length)
                pos[j].x = length - body_radius;
            
            if ((pos[j].y - body_radius) < 0)
                pos[j].y = body_radius;
            else if ((pos[j].y + body_radius) > width)
                pos[j].y = width - body_radius;

            if ((pos[j].z - body_radius) < 0)
                pos[j].z = body_radius;
            else if ((pos[j].z + body_radius) > depth)
                pos[j].z = depth - body_radius;
            
            vel[j].x = vhs.x + (force[j].x*time_step)/(2*body_mass);
            vel[j].y = vhs.y+ (force[j].y*time_step)/(2*body_mass);
            vel[j].z = vhs.z + (force[j].z*time_step)/(2*body_mass);

            //collision with the wall
            if ((pos[j].x - body_radius) <= 0 || (pos[j].x + body_radius) >= length){
                vel[j].x = -vel[j].x;
                
            }

            if ((pos[j].y - body_radius) <= 0 || (pos[j].y + body_radius) >= width){
                vel[j].y = -vel[j].y;   
            }

            if ((pos[j].z - body_radius) <= 0 || (pos[j].z + body_radius) >= depth){
                vel[j].z = -vel[j].z;
            }

        }
        #pragma omp barrier

        #pragma omp single
        for(j=0; j<num_bodies; j++){
            if ((i+1)%100 == 0){
                fprintf(ptr, "%0.2lf %0.2lf %0.2lf\n", pos[j].x, pos[j].y, pos[j].z);
            }
        } 
        #pragma omp for 
        for(j=0; j<num_bodies - 1; j++){
            for(k=j+1; k<num_bodies; k++){
                r.x = pos[k].x - pos[j].x;
                r.y = pos[k].y - pos[j].y;
                r.z = pos[k].z - pos[j].z;
                dist = sqrt(r.x*r.x + r.y*r.y + r.z*r.z);
                //collision with other particles
                if (dist <= 2*body_radius){
                    temp = vel[j];
                    vel[j] = vel[k];
                    vel[k] = temp;
                }
            }
        }
        #pragma omp barrier
    }
    if((i+1)%1000 == 0)
    printf("Iteration %d completed!\n", i+1);
    }
    wtime = omp_get_wtime() - wtime;
    printf("Time elapsed: %lf", wtime);
}

int main(){
    FILE *ptr;
    char ch;
    ptr = fopen("Trajectory.txt", "r");
    double length, width, depth, body_radius, body_mass, time_step;
    int num_bodies;
    fscanf(ptr, "%*s %*s %*s %*s %*s %lf", &length);
    fscanf(ptr, "%*s %*s %*s %*s %*s %lf", &width);
    fscanf(ptr, "%*s %*s %*s %*s %*s %lf", &depth);
    fscanf(ptr, "%*s %*s %*s %*s %*s %d", &num_bodies);
    fscanf(ptr, "%*s %*s %*s %*s %*s %lf", &body_radius);
    fscanf(ptr, "%*s %*s %*s %*s %*s %lf", &body_mass);
    fscanf(ptr, "%*s %*s %*s %*s %*s %*s %*s %*s %*s %lf", &time_step);

    fscanf(ptr, "%*s %*s %*s");
    
    vector pos[1000], vel[1000], force[1000];
    for(int i=0; i<1000; i++){
        fscanf(ptr, "%lf %lf %lf", &pos[i].x, &pos[i].y, &pos[i].z);
        vel[i].x = 0;
        vel[i].y = 0;
        vel[i].z = 0;
        force[i].x = 0;
        force[i].y = 0;
        force[i].z = 0;
    }
    int n_sims = 720000;

    get_trajectory(length, width, depth, num_bodies, body_radius, body_mass, time_step, pos, vel, force, n_sims);
}