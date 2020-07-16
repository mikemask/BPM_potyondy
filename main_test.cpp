#include <iomanip>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include <list>
#include "abs.h"
#include "part.h"
#include "distance.h"
#include "prop_bond_potyondy.h"
#include "forces_bond_potyondy.h"

int main()
{
    double const rho = 7000.;
    double const d = 0.001;
    double const E = 2.7e6;
    double const nu = 0.4;
    double const dt = 0.0000001;
    double const eps = d/1000.;
    double const lambda = 1.1;
    double const gamma = 1.;
    int bonds = 0;
    int count = 0;
    int const n = 10000000;
    int const n_part = 2;

    std::ofstream out_file;
    out_file.open("data_test_hyb_velx.csv");
    out_file << "timestep,ID,x,y,z,u,v,w,fx,fy,fz,omx,omy,omz,tx,ty,tz\n";
    int const n_out = 10000;
    bool flag = false;

    Part vec_part[n_part];

    Bond vec_bond;

    /*Initial particle configuration*/

    double const m = rho*(4./3.)*3.14159265*pow((d/2.),3.);

    double pos1[3], pos2[3];

    pos1[0] = 0.;
    pos1[1] = 0.;
    pos1[2] = 0.;

    pos2[0] = 0.;
    pos2[1] = 0.;
    pos2[2] = d+eps;

    vec_part[0].setPos(pos1);
    vec_part[1].setPos(pos2);
    vec_part[0].setPosOld(pos1);
    vec_part[1].setPosOld(pos2);

    for (int i=0; i<n_part; i++)
    {
        vec_part[i].setMass(m);
        vec_part[i].setDiameter(d);
    }


//    for (int i=0; i<n_part; i++)
//    {
//        std::cout << "particle:" << i << "  " << "X:" << vec_part[i].getPos()[0] << ","
//               << vec_part[i].getPos()[1] << "," << vec_part[i].getPos()[2] <<"\n" ;
//    }


    /*Bond Initialization*/

    double *posi, *posj, radi, radj, radi_cr, radj_cr, dist;
    int IDS[2];

    posi = vec_part[0].getPos();
    radi = vec_part[0].getDiameter()/2.;
    radi_cr = radi*lambda;

    posj = vec_part[1].getPos();
    radj = vec_part[1].getDiameter()/2.;
    radj_cr = radj*lambda;

    dist = distance(posi, posj);

    if (dist < (radi_cr+radj_cr))
    {
        double radius;
        double const kn = E/(radi+radj);
        double const ks = kn/2.5;
        double inertia[2];
//        double const kn = sqrt(2.)*E*radi/(2.*(1.-2.*nu));
//        double const ks = kn/2.5;

        radius = lambda*std::min(radi,radj);

        IDS[0] = 0;
        IDS[1] = 1;

        inertia[0] = 0.25*3.14159265*pow(radius,4.);
        inertia[1] = 2.*inertia[0];

        vec_bond.setRad(radius);
        vec_bond.setStiffNorm(kn);
        vec_bond.setStiffShear(ks);
        vec_bond.setIds(IDS);
        vec_bond.setInertia(inertia);
    }


    /*Starts the dynamics*/

    for (int i=0; i<n; i++)
    {
        /*first timestep, initialization*/

        if(i==0)
        {
            double *x, *v, *om;

            x = vec_part[1].getPos();
            v = vec_part[1].getVel();
            om = vec_part[1].getRot();

//            x[0] += 1.e-6;
//            x[2] += 1.e-6;
            v[0] = 0.001;
//            v[2] = 0.001;
//            om[0] = -0.0001;
//            om[1] = 0.001;
        }

        else
        {
            /*New positions with integration of EOM due to bond forces*/

            for (int j=0; j<n_part; j++)
            {
                double *x, *v, *f, *om, *t;
                double *v_old, *f_old, *t_old;
                double x_new[3], v_new[3], om_new[3];
                double x_act[3], v_act[3];
                double m, d, inertia;

                x = vec_part[j].getPos();
                v = vec_part[j].getVel();
                f = vec_part[j].getForce();
                om = vec_part[j].getRot();
                t = vec_part[j].getTorque();

                v_old = vec_part[j].getVelOld();
                f_old = vec_part[j].getForceOld();
                t_old = vec_part[j].getTorqueOld();

                m = vec_part[j].getMass();
                d = vec_part[j].getDiameter();
                inertia = (2./5.)*m*pow((d/2.),2.);

                if(i==1)
                {
                    for (int k=0; k<3; k++)
                    {
                        x_act[k] = x[k];  
                        v_act[k] = v[k];

                        x_new[k] = x[k] + v[k]*dt;
                        v_new[k] = v[k] + (dt/m)*f[k];
                        om_new[k] = om[k] + (dt/inertia)*t[k];
                    }
                }

                else
                {
                    for (int k=0; k<3; k++)
                    {
                        x_act[k] = x[k];
                        v_act[k] = v[k];

                        x_new[k] = x[k] + 0.5*dt*(3.*v[k] - v_old[k]);
                        v_new[k] = v[k] + 0.5*(dt/m)*(3.*f[k] - f_old[k]);
                        om_new[k] = om[k] + 0.5*(dt/inertia)*(3.*t[k] - t_old[k]);
                    }
                }

                vec_part[j].setPos(x_new);
                vec_part[j].setVel(v_new);
                vec_part[j].setRot(om_new);

                vec_part[j].setPosOld(x_act);
                vec_part[j].setVelOld(v_act);
                vec_part[j].setForceOld(f);
                vec_part[j].setTorqueOld(t);
            }
        }

        //Forces on bonds

        Bond *bondino;

        bondino = &vec_bond;
        forces_bond(bondino, vec_part[0], vec_part[1], dt, gamma);

        if (i==0)
            std::cout << "normal force: " << vec_bond.getForce_n()[0] << "," << vec_bond.getForce_n()[1] << "," << vec_bond.getForce_n()[2] << "\n"
                      << "shear-trans force: " << vec_bond.getForce_s()[0] << "," << vec_bond.getForce_s()[1] << "," << vec_bond.getForce_s()[2] << "\n"
                      << "shear-trans torque: " << vec_bond.getTorque_st()[0] << "," << vec_bond.getTorque_st()[1] << "," << vec_bond.getTorque_st()[2] << "\n";

        //Forces on particles

        for (int j=0; j<n_part; j++)
        {
            double F_bond[3], T_bond[3];

            //Bonds contribution

            int *IDS;

            IDS = vec_bond.getIds();

            if (j == IDS[0])
            {
                F_bond[0] = -(vec_bond.getForce_n()[0] + vec_bond.getForce_s()[0]);
                F_bond[1] = -(vec_bond.getForce_n()[1] + vec_bond.getForce_s()[1]);
                F_bond[2] = -(vec_bond.getForce_n()[2] + vec_bond.getForce_s()[2]);

                T_bond[0] = -(vec_bond.getTorque_s()[0] + vec_bond.getTorque_st()[0] + vec_bond.getTorque_n()[0]);
                T_bond[1] = -(vec_bond.getTorque_s()[1] + vec_bond.getTorque_st()[1] + vec_bond.getTorque_n()[1]);
                T_bond[2] = -(vec_bond.getTorque_s()[2] + vec_bond.getTorque_st()[2] + vec_bond.getTorque_n()[2]);
            }

            else if (j == IDS[1])
            {
                F_bond[0] = (vec_bond.getForce_n()[0] + vec_bond.getForce_s()[0]);
                F_bond[1] = (vec_bond.getForce_n()[1] + vec_bond.getForce_s()[1]);
                F_bond[2] = (vec_bond.getForce_n()[2] + vec_bond.getForce_s()[2]);

                T_bond[0] = vec_bond.getTorque_s()[0] + vec_bond.getTorque_st()[0] + vec_bond.getTorque_n()[0];
                T_bond[1] = vec_bond.getTorque_s()[1] + vec_bond.getTorque_st()[1] + vec_bond.getTorque_n()[1];
                T_bond[2] = vec_bond.getTorque_s()[2] + vec_bond.getTorque_st()[2] + vec_bond.getTorque_n()[2];
            }


            vec_part[j].setForce(F_bond);
            vec_part[j].setTorque(T_bond);
        }

        for(int j=0; j<n_part; j++)
        {
            if (i == 10 || i == 20)
            {
                std::cout << j << ","
                          << vec_part[j].getPos()[0] << "," << vec_part[j].getPos()[1] << "," << vec_part[j].getPos()[2]
                          << "," << vec_part[j].getVel()[0] << "," << vec_part[j].getVel()[1] << "," << vec_part[j].getVel()[2]
                          << "," << vec_part[j].getForce()[0] << "," << vec_part[j].getForce()[1] << "," << vec_part[j].getForce()[2]
                          << "," << vec_part[j].getTorque()[0] << "," << vec_part[j].getTorque()[1] << "," << vec_part[j].getTorque()[2] << "\n";
            }
        }

        for(int j=0; j<n_part; j++)
        {
            if (i == count*n_out)
            {
                out_file << i*dt << "," << j << ","
                         << vec_part[j].getPos()[0] << "," << vec_part[j].getPos()[1] << "," << vec_part[j].getPos()[2] << ","
                         << vec_part[j].getVel()[0] << "," << vec_part[j].getVel()[1] << "," << vec_part[j].getVel()[2] << ","
                         << vec_part[j].getForce()[0] << "," << vec_part[j].getForce()[1] << "," << vec_part[j].getForce()[2] << ","
                         << vec_part[j].getRot()[0] << "," << vec_part[j].getRot()[1] << "," << vec_part[j].getRot()[2] << ","
                         << vec_part[j].getTorque()[0] << "," << vec_part[j].getTorque()[1] << "," << vec_part[j].getTorque()[2] << "\n";
                flag = true;
            }
        }

        if (flag)
        {
            count++;
            flag = false;
        }
    }

    for (int i=0; i<n_part; i++)
    {
        std::cout << i << "," << vec_part[i].getPos()[0] << ","
                              << vec_part[i].getPos()[1] << ","
                              << vec_part[i].getPos()[2] << "\n";
    }

    return 0;
}
