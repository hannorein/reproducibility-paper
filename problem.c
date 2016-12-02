#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <sys/stat.h>
#include "rebound.h"


void gr_force(struct reb_simulation* const r){
    // From REBOUNDx
    const struct reb_particle source = r->particles[0];
    const double C2 = 1.0130251e+08;
    const double prefac1 = 6.*(r->G*source.m)*(r->G*source.m)/C2;
    
	for (int i=1;i<r->N;i++){
        const struct reb_particle p = r->particles[i];
        const double dx = p.x - source.x;
        const double dy = p.y - source.y;
        const double dz = p.z - source.z;
        const double r2 = dx*dx + dy*dy + dz*dz;
        const double prefac = prefac1/(r2*r2);
        r->particles[i].ax -= prefac*dx;
        r->particles[i].ay -= prefac*dy;
        r->particles[i].az -= prefac*dz;
        r->particles[0].ax += p.m/source.m*prefac*dx;
        r->particles[0].ay += p.m/source.m*prefac*dy;
        r->particles[0].az += p.m/source.m*prefac*dz;
    }
}

int main(int argc, char* argv[]) {
    int procs_N = 1;
    int procs_N_start = 0;
    if (argc>=2){
        procs_N = atoi(argv[1]);    
    }
    if (argc>=3){
        procs_N_start = atoi(argv[2]);    
    }
    int proc_i;
    for(proc_i=procs_N_start; proc_i<procs_N-1; proc_i++) {
        int pid = fork();
        if (pid == 0) {
            break;
        }
    }

    // Read initial conditions

    char filename[512];
    sprintf(filename,"restart_%04d.bin",proc_i);
    struct reb_simulation* r = reb_simulationarchive_restart(filename);
    if (r==NULL){
        r= reb_create_simulation_from_binary("../ss-2016-06-18.bin");
        r->particles[1].x += 6.6845871e-12*(double)proc_i; // Shift by 1 m * proc_i
        reb_move_to_com(r);
        r->dt = 6./365.25*2.*M_PI;    
        r->simulationarchive_interval = 2.*M_PI*5.1434563422923e4;
        r->ri_whfast.safe_mode = 0;     
        r->ri_whfast.corrector = 5;    
        r->integrator = REB_INTEGRATOR_WHFAST;
        r->gravity = REB_GRAVITY_BASIC;
    }else{
        printf("Loaded simulation. Time now: %.16f\n",r->t);
    }

    // Setup simulation
    r->simulationarchive_filename = filename;
    r->additional_forces = gr_force;

    reb_integrate(r, 2.*M_PI*1e10); //10 Gyr
}
