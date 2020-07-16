#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include <list>
#include "abs.h"
#include "part.h"
#include "distance.h"
#include "forces_bond_potyondy.h"
#include "prop_bond_potyondy.h"

int main()
{
    double const rho = 7000.;
    double const d = 0.001;
    double const E = 2.7e9;
    double const nu = 0.4;
    double const dt = 0.0000001;
    double const T = 0.1;
    double const eps = d/1000.;
    double const lambda = 1.1;
    double const gamma = 1.;
    int bonds = 0;
    int count = 0;
    int const n = 10000000;
    int const n_part = 27;

    std::ofstream out_file;
    out_file.open("dataset.csv");
    out_file << "timestep,ID,x,y,z,u,v,w,omx,omy,omz,fx,fy,fz,tx,ty,tz \n";
    int n_out = 10000;
    bool flag = false;

    Part vec_part[n_part];

    Bond *vec_bond = new Bond[100];

////Initial particle configuration

    double const m = rho*(4./3.)*3.14159265*pow((d/2.),3.);

    int id = 0;
    double pos[3]={};

    for (int i=0; i<3; i++)
    {
        for (int j=0; j<3; j++)
        {
            for (int k=0; k<3; k++)
            {
                pos[0] = k*(d+eps);
                pos[1] = j*(d+eps);
                pos[2] = i*(d+eps);

                vec_part[id].setMass(m);
                vec_part[id].setDiameter(d);
                vec_part[id].setPos(pos);
                vec_part[id].setPosIn(pos);
                id++;
            }
        }
    }

    for (int i=0; i<n_part; i++)
    {
        std::cout << "particle:" << i << "  " << "X:" << vec_part[i].getPos()[0] << ","
               << vec_part[i].getPos()[1] << "," << vec_part[i].getPos()[2] <<"\n" ;
    }

////Bond initialization

    int pp = 1;

    for (int i=0; i<n_part; i++)
    {
        for (int j=pp; j<n_part; j++)
        {
            if (i != j)
            {
                double *posi, *posj, radi, radj, radi_cr, radj_cr, dist;
                int IDS[2];

                posi = vec_part[i].getPos();
                radi = vec_part[i].getDiameter()/2.;
                radi_cr = radi*lambda;

                posj = vec_part[j].getPos();
                radj = vec_part[j].getDiameter()/2.;
                radj_cr = radj*lambda;

                dist = distance(posi, posj);

                if (dist <= (radi_cr+radj_cr))
                {
                    double radius;
                    //double const kn = E/(radi+radj);
                    //double const ks = kn/2.5;
                    double const kn = sqrt(2.)*E*radi/(2.*(1.-2.*nu));
                    double const ks = kn*(1.-3.*nu)/(1.+nu);

                    radius = lambda*std::min(radi,radj);

                    IDS[0] = i;
                    IDS[1] = j;

                    vec_bond[bonds].setRad(radius);
                    vec_bond[bonds].setStiffNorm(kn);
                    vec_bond[bonds].setStiffShear(ks);
                    vec_bond[bonds].setIds(IDS);

                    bonds++;
                }
            }
        }

        pp++;
    }

    /*Starts the dynamics*/

    for (int i=0; i<n; i++)
    {
        /*first timestep, initialization*/

        if(i==0)
        {
            for (int j=18; j<n_part; j++)
            {
                double *pos;

                pos = vec_part[j].getPos();

                pos[0] += 1.e-8;

                vec_part[j].setPos(pos);

                //Forces on bonds

                for (int j=0; j<bonds; j++)
                {
                    int *IDS;
                    Bond *bondino;

                    IDS = vec_bond[j].getIds();
                    bondino = &vec_bond[j];
                    bond_fdm(bondino, vec_part[IDS[0]], vec_part[IDS[1]], dt);
//                    f_bonds(bondino, vec_part[IDS[0]], vec_part[IDS[1]], dt, gamma);
                }
            }

            //Forces on particles

            for (int j=0; j<n_part; j++)
            {
                double F[3]={}, T[3]={};
                double F_bond[3]={}, T_bond[3]={};

                //Bonds contribution

                for (int k=0; k<bonds; k++)
                {
                    int *IDS;

                    IDS = vec_bond[k].getIds();

                    double q_conj[4]={};
                    double *q_f;

                    if (j == IDS[0])
                    {
                        double *q = vec_part[IDS[1]].getQuat();

                        F_bond[0] = vec_bond[k].getForce_n()[0] + vec_bond[k].getForce_st()[0];
                        F_bond[1] = vec_bond[k].getForce_n()[1] + vec_bond[k].getForce_st()[1];
                        F_bond[2] = vec_bond[k].getForce_n()[2] + vec_bond[k].getForce_st()[2];

                        //q_conj[0] = q[0];
                        //q_conj[1] = -q[1];
                        //q_conj[2] = -q[2];
                        //q_conj[3] = -q[3];

                        //q_f = rotation(q_conj, F_bond);

                    }

                    else if (j == IDS[1])
                    {
                        double *q = vec_part[IDS[0]].getQuat(); //check this might be wrong

                        F_bond[0] = -vec_bond[k].getForce_n()[0] - vec_bond[k].getForce_st()[0];
                        F_bond[1] = -vec_bond[k].getForce_n()[1] - vec_bond[k].getForce_st()[1];
                        F_bond[2] = -vec_bond[k].getForce_n()[2] - vec_bond[k].getForce_st()[2];

                        //q_conj[0] = q[0];
                        //q_conj[1] = -q[1];
                        //q_conj[2] = -q[2];
                        //q_conj[3] = -q[3];

                        //q_f = rotation(q_conj, F_bond);

                    }

                    F[0] += F_bond[0];
                    F[1] += F_bond[1];
                    F[2] += F_bond[2];
                }
                //New positions and velocities from EOM with Forward Euler

                double *x, *v, *om, *f, *t;
                double v_act[3];
                double m, d, inertia;

                x = vec_part[j].getPos();
                v = vec_part[j].getVel();
                om = vec_part[j].getRot();
                f = vec_part[j].getForce();

                m = vec_part[j].getMass();
                d = vec_part[j].getDiameter();
                inertia = (2./5.)*m*pow((d/2.),2.);

                for (int k=0; k<3; k++)
                {
                    v_act[k] = v[k];

                    x[k] += v[k]*dt;
                    v[k] += f[k]*dt/m;
                }
            }
        }

        else
        {
            //Forces on bonds

            for (int j=0; j<bonds; j++)
            {
                int *IDS;
                Bond *bondino;

                IDS = vec_bond[j].getIds();
                bondino = &vec_bond[j];
                bond_fdm(bondino, vec_part[IDS[0]], vec_part[IDS[1]], dt);
//                f_bonds(bondino, vec_part[IDS[0]], vec_part[IDS[1]], dt, gamma);
            }

            //Forces on particles

            for (int j=0; j<n_part; j++)
            {
                double F[3]={}, T[3]={};
                double F_bond[3]={}, T_bond[3]={};

                //Bonds contribution

                for (int k=0; k<bonds; k++)
                {
                    int *IDS;

                    IDS = vec_bond[k].getIds();

                    double q_conj[4]={};
                    double *q_f;

                    if (j == IDS[0])
                    {
                        double *q = vec_part[IDS[1]].getQuat();

                        F_bond[0] = vec_bond[k].getForce_n()[0] + vec_bond[k].getForce_st()[0];
                        F_bond[1] = vec_bond[k].getForce_n()[1] + vec_bond[k].getForce_st()[1];
                        F_bond[2] = vec_bond[k].getForce_n()[2] + vec_bond[k].getForce_st()[2];

                        //q_conj[0] = q[0];
                        //q_conj[1] = -q[1];
                        //q_conj[2] = -q[2];
                        //q_conj[3] = -q[3];

                        //q_f = rotation(q_conj, F_bond);

                    }

                    else if (j == IDS[1])
                    {
                        double *q = vec_part[IDS[0]].getQuat();

                        F_bond[0] = -vec_bond[k].getForce_n()[0] - vec_bond[k].getForce_st()[0];
                        F_bond[1] = -vec_bond[k].getForce_n()[1] - vec_bond[k].getForce_st()[1];
                        F_bond[2] = -vec_bond[k].getForce_n()[2] - vec_bond[k].getForce_st()[2];

                        //q_conj[0] = q[0];
                        //q_conj[1] = -q[1];
                        //q_conj[2] = -q[2];
                        //q_conj[3] = -q[3];

                        //q_f = rotation(q_conj, F_bond);

                    }

                    F[0] += F_bond[0];
                    F[1] += F_bond[1];
                    F[2] += F_bond[2];
                }

                vec_part[j].setForce(F);

                //New positions and velocities from EOM with Adams-Bashforth

                double *x, *v, *om, *f, *t, *v_old, *f_old, *t_old;
                double v_act[3];
                double m, d, inertia;

                x = vec_part[j].getPos();
                v = vec_part[j].getVel();
                om = vec_part[j].getRot();
                f = vec_part[j].getForce();

                v_old = vec_part[j].getVelOld();
                f_old = vec_part[j].getForceOld();

                m = vec_part[j].getMass();
                d = vec_part[j].getDiameter();
                inertia = (2./5.)*m*pow((d/2.),2.);

                for (int k=0; k<3; k++)
                {
                    v_act[k] = v[k];

                    x[k] += 0.5*dt*(3.*v[k] - v_old[k]);
                    v[k] += 0.5*(dt/m)*(3.*f[k] - f_old[k]);
                }

                vec_part[j].setPos(x);
                vec_part[j].setVel(v);

                vec_part[j].setVelOld(v_act);
                vec_part[j].setForceOld(f);

                if (i==1 || i==2)
                {
                    std::cout << j << "," << vec_part[j].getForce()[0] << "," << vec_part[j].getForce()[1] << "," << vec_part[j].getForce()[2] << "\n";
                }

                T_0 = i*dt;

                if (fabs(T_0-count*T_out) < 1e-5)
                {
                    out_file << i*dt << "," << j << ","
                              << vec_part[j].getPos()[0] << "," << vec_part[j].getPos()[0] << "," << vec_part[j].getPos()[2]
                              << "," << vec_part[j].getVel()[0] << "," << vec_part[j].getVel()[1] << "," << vec_part[j].getVel()[2]
                              << "," << vec_part[j].getRot()[0] << "," << vec_part[j].getRot()[1] << "," << vec_part[j].getRot()[2]
                              << "," << vec_part[j].getForce()[0] << "," << vec_part[j].getForce()[1] << "," << vec_part[j].getForce()[2]
                              << "," << vec_part[j].getTorque()[0] << "," << vec_part[j].getTorque()[1] << "," << vec_part[j].getTorque()[2] << "\n";
                    flag = true;
                }
            }

            if (flag)
            {
                count++;
                flag = false;
                std::cout << count << "\n";
            }
        }
    }

    for (int i=0; i<n_part; i++)
    {
        std::cout << i << "," << vec_part[i].getPos()[0] << ","
                              << vec_part[i].getPos()[1] << ","
                              << vec_part[i].getPos()[2] << "\n";
    }

    delete [] vec_bond;
    return 0;
}
