void forces_bond(Bond *bond, Part parti, Part partj, double dt, double gamma)
{

//Getting bonds values

    double kn, ks, *fn, *fs, *tn, *ts, *tst, *inertia, radius, area;    

    kn = bond -> getStiffNorm();
    ks = bond -> getStiffShear();

    fn = bond -> getForce_n();
    fs = bond -> getForce_s();
    tn = bond -> getTorque_n();
    ts = bond -> getTorque_s();
    tst = bond -> getTorque_st();

    radius = bond -> getRad();
    area = 3.14159265*radius*radius;
    inertia = bond -> getInertia();

//Getting particles values

    double *vel_i, *vel_j, *om_i, *om_j, *pos_i, *pos_j, *pos_i_old, *pos_j_old;

    vel_i = parti.getVel();
    vel_j = partj.getVel();

    om_i = parti.getRot();
    om_j = partj.getRot();

    pos_i = parti.getPos();
    pos_j = partj.getPos();

    pos_i_old = parti.getPosOld();
    pos_j_old = partj.getPosOld();

//Relative motion

    double v_cont[3];
    double v_rel[3];
    double om_rel[3];
    double x_cont[3];
    double N[3], n[3], abs_N, N_old[3], n_old[3], abs_Nold;

    v_rel[0] = vel_i[0] - vel_j[0];
    v_rel[1] = vel_i[1] - vel_j[1];
    v_rel[2] = vel_i[2] - vel_j[2];

    om_rel[0] = om_i[0] - om_j[0];
    om_rel[1] = om_i[1] - om_j[1];
    om_rel[2] = om_i[2] - om_j[2];

    x_cont[0] = 0.5*(pos_i[0] + pos_j[0]);
    x_cont[1] = 0.5*(pos_i[1] + pos_j[1]);
    x_cont[2] = 0.5*(pos_i[2] + pos_j[2]);

    N[0] = pos_i[0] - pos_j[0];
    N[1] = pos_i[1] - pos_j[1];
    N[2] = pos_i[2] - pos_j[2];

    N_old[0] = pos_i_old[0] - pos_j_old[0];
    N_old[1] = pos_i_old[1] - pos_j_old[1];
    N_old[2] = pos_i_old[2] - pos_j_old[2];
    
    abs_N = abs(N);
    abs_Nold = abs(N_old);

    n[0] = N[0]/abs_N;
    n[1] = N[1]/abs_N;
    n[2] = N[2]/abs_N;

    n_old[0] = N_old[0]/abs_Nold;
    n_old[1] = N_old[1]/abs_Nold;
    n_old[2] = N_old[2]/abs_Nold;

    double t[3];

    t[0] = n_old[1]*n[2] - n_old[2]*n[1];
    t[1] = n_old[2]*n[0] - n_old[0]*n[2];
    t[2] = n_old[0]*n[1] - n_old[1]*n[0];

    double s[3];

    s[0] = t[1]*n[2] - t[2]*n[1];
    s[1] = t[2]*n[0] - t[0]*n[2];
    s[2] = t[0]*n[1] - t[1]*n[0];
    
    double alfa;

    alfa = acos(n_old[0]*n[0] + n_old[1]*n[1] + n_old[2]*n[2]);

    int e;

    for (int i=0; i<3; i++)
    {
        double perm = 0.;

        for (int j=0; j<3; j++)
        {
            for (int k=0; k<3; k++)
            {
                if (i==j || i==k || j==k || i==j==k)
                    e = 0;
                else if ((i==1 && j==2 && k==3) || (i==2 && j==3 && k==1) || (i==3 && j==1 && k==2))
                    e = 1;
                else
                    e = -1;

                perm += e*(om_i[j]*(x_cont[k] - pos_i[k]) - om_j[j]*(x_cont[k] - pos_j[k]));
            }
        }

        v_cont[i] = v_rel[i] + perm;
    }

    double v_cont_n[3]={}, v_cont_s[3]={};
    double om_rel_n[3]={}, om_rel_s[3]={};
    double delta_un[3], delta_us[3];
    double delta_theta_n[3], delta_theta_s[3], abs_thetan;

    for (int i=0; i<3; i++)
    {
        for (int j=0; j<3; j++)
        {
            v_cont_n[i] += v_cont[j]*n[i]*n[j];
            om_rel_n[i] += om_rel[j]*n[i]*n[j];

            v_cont_s[i] += v_cont[j]*s[i]*s[j];
            om_rel_s[i] += om_rel[j]*s[i]*s[j];
        }
    }

    delta_un[0] = v_cont_n[0]*dt;
    delta_un[1] = v_cont_n[1]*dt;
    delta_un[2] = v_cont_n[2]*dt;

    delta_us[0] = v_cont_s[0]*dt;
    delta_us[1] = v_cont_s[1]*dt;
    delta_us[2] = v_cont_s[2]*dt;

    delta_theta_n[0] = om_rel_n[0]*dt;
    delta_theta_n[1] = om_rel_n[1]*dt;
    delta_theta_n[2] = om_rel_n[2]*dt;

    delta_theta_s[0] = om_rel_s[0]*dt;
    delta_theta_s[1] = om_rel_s[1]*dt;
    delta_theta_s[2] = om_rel_s[2]*dt;

    abs_thetan = abs(delta_theta_n);

//Force increment computation and total force update

    double delta_fn[3], delta_fs[3], abs_fs, fs_rot[3], delta_tst[3];

    delta_fn[0] = kn*area*delta_un[0];
    delta_fn[1] = kn*area*delta_un[1];
    delta_fn[2] = kn*area*delta_un[2];

    delta_fs[0] = ks*area*delta_us[0];
    delta_fs[1] = ks*area*delta_us[1];
    delta_fs[2] = ks*area*delta_us[2];

    abs_fs = abs(delta_fs);

    fn[0] = fn[0]*gamma + delta_fn[0];
    fn[1] = fn[1]*gamma + delta_fn[1];
    fn[2] = fn[2]*gamma + delta_fn[2];

    fs_rot[0] = fs[0]*cos(abs_thetan) + (n[1]*fs[2] - n[2]*fs[1])*sin(abs_thetan) + n[0]*(n[0]*fs[0] + n[1]*fs[1] + n[2]*fs[2])*(1-cos(abs_thetan));
    fs_rot[1] = fs[1]*cos(abs_thetan) + (n[2]*fs[0] - n[0]*fs[2])*sin(abs_thetan) + n[1]*(n[0]*fs[0] + n[1]*fs[1] + n[2]*fs[2])*(1-cos(abs_thetan));
    fs_rot[2] = fs[2]*cos(abs_thetan) + (n[0]*fs[1] - n[1]*fs[0])*sin(abs_thetan) + n[2]*(n[0]*fs[0] + n[1]*fs[1] + n[2]*fs[2])*(1-cos(abs_thetan));

    fs[0] = fs_rot[0]*gamma + delta_fs[0];
    fs[1] = fs_rot[1]*gamma + delta_fs[1];
    fs[2] = fs_rot[2]*gamma + delta_fs[2];

    delta_tst[0] = 0.5*abs_N*abs_fs*t[0];
    delta_tst[1] = 0.5*abs_N*abs_fs*t[1];
    delta_tst[2] = 0.5*abs_N*abs_fs*t[2];

    tst[0] += delta_tst[0];
    tst[1] += delta_tst[1];
    tst[2] += delta_tst[2];

}
