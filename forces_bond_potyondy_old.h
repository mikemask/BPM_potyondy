void f_bonds(Bond *bond, Part parti, Part partj, double dt, double gamma)
{

//Getting bonds values

    double kn, ks, *fn, *fs, *tn, *ts, *tst, *inertia;    

    kn = bond -> getStiffNorm();
    ks = bond -> getStiffShear();

    fn = bond -> getForce_n();
    fs = bond -> getForce_st();
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

//Normal versor computation

    norm_mod = distance(pos_i,pos_j);

    norm[0] = pos_j[0]-pos_i[0];
    norm[1] = pos_j[1]-pos_i[1];
    norm[2] = pos_j[2]-pos_i[2];

    norm_vers[0] = norm[0]/norm_mod;
    norm_vers[1] = norm[1]/norm_mod;
    norm_vers[2] = norm[2]/norm_mod;

//Vectors projections into normal and shear direction w.r.t. contact plane

    scal_vel = rel_vel[0]*norm_vers[0]+rel_vel[1]*norm_vers[1]+rel_vel[2]*norm_vers[2];
    scal_rot = rel_rot[0]*norm_vers[0]+rel_rot[1]*norm_vers[1]+rel_rot[2]*norm_vers[2];

    rel_vel_n[0] = scal_vel*norm_vers[0];
    rel_vel_n[1] = scal_vel*norm_vers[1];
    rel_vel_n[2] = scal_vel*norm_vers[2];

    rel_rot_n[0] = scal_rot*norm_vers[0];
    rel_rot_n[1] = scal_rot*norm_vers[1];
    rel_rot_n[2] = scal_rot*norm_vers[2];

//    rel_vel_s[0] = rel_vel[0] - rel_vel_n[0];
//    rel_vel_s[1] = rel_vel[1] - rel_vel_n[1];
//    rel_vel_s[2] = rel_vel[2] - rel_vel_n[2];
//
//    rel_rot_s[0] = rel_rot[0] - rel_rot_n[0];
//    rel_rot_s[1] = rel_rot[1] - rel_rot_n[1];
//    rel_rot_s[2] = rel_rot[2] - rel_rot_n[2];

//Computation of rotation angle and rotational matrix

    theta = sqrt((rel_rot_n[0]*rel_rot_n[0])+(rel_rot_n[1]*rel_rot_n[1])+(rel_rot_n[2]*rel_rot_n[2]))*dt;

    R[0][0] = cos(theta) + norm_vers[0]*norm_vers[0]*(1-cos(theta));
    R[0][1] = norm_vers[0]*norm_vers[1]*(1-cos(theta))-norm_vers[2]*sin(theta);
    R[0][2] = norm_vers[0]*norm_vers[2]*(1-cos(theta))+norm_vers[1]*sin(theta);

    R[1][0] = norm_vers[0]*norm_vers[1]*(1-cos(theta))+norm_vers[2]*sin(theta);
    R[1][1] = cos(theta) + norm_vers[1]*norm_vers[1]*(1-cos(theta));
    R[1][2] = norm_vers[1]*norm_vers[2]*(1-cos(theta))-norm_vers[0]*sin(theta);

    R[2][0] = norm_vers[0]*norm_vers[2]*(1-cos(theta))-norm_vers[1]*sin(theta);
    R[2][1] = norm_vers[1]*norm_vers[2]*(1-cos(theta))+norm_vers[0]*sin(theta);
    R[2][2] = cos(theta) + norm_vers[2]*norm_vers[2]*(1-cos(theta));

//Computation of rotated shear vectors

    rel_vel_s[0] = R[0][0]*(rel_vel[0]-rel_vel_n[0])+R[0][1]*(rel_vel[1]-rel_vel_n[1])+R[0][2]*(rel_vel[2]-rel_vel_n[2]);
    rel_vel_s[1] = R[1][0]*(rel_vel[0]-rel_vel_n[0])+R[1][1]*(rel_vel[1]-rel_vel_n[1])+R[1][2]*(rel_vel[2]-rel_vel_n[2]);
    rel_vel_s[2] = R[2][0]*(rel_vel[0]-rel_vel_n[0])+R[2][1]*(rel_vel[1]-rel_vel_n[1])+R[2][2]*(rel_vel[2]-rel_vel_n[2]);

    rel_rot_s[0] = R[0][0]*(rel_rot[0]-rel_rot_n[0])+R[0][1]*(rel_rot[1]-rel_rot_n[1])+R[0][2]*(rel_rot[2]-rel_rot_n[2]);
    rel_rot_s[1] = R[1][0]*(rel_rot[0]-rel_rot_n[0])+R[1][1]*(rel_rot[1]-rel_rot_n[1])+R[1][2]*(rel_rot[2]-rel_rot_n[2]);
    rel_rot_s[2] = R[2][0]*(rel_rot[0]-rel_rot_n[0])+R[2][1]*(rel_rot[1]-rel_rot_n[1])+R[2][2]*(rel_rot[2]-rel_rot_n[2]);

//Force increment computation and total force update

    delta_fn[0] = kn*area*dt*rel_vel_n[0];
    delta_fn[1] = kn*area*dt*rel_vel_n[1];
    delta_fn[2] = kn*area*dt*rel_vel_n[2];

    delta_fs[0] = -ks*area*dt*rel_vel_s[0];
    delta_fs[1] = -ks*area*dt*rel_vel_s[1];
    delta_fs[2] = -ks*area*dt*rel_vel_s[2];

    fn[0] = fn[0]*gamma + delta_fn[0];
    fn[1] = fn[1]*gamma + delta_fn[1];
    fn[2] = fn[2]*gamma + delta_fn[2];

    fs[0] = fs[0]*gamma + delta_fs[0];
    fs[1] = fs[1]*gamma + delta_fs[1];
    fs[2] = fs[2]*gamma + delta_fs[2];


    bond -> setForce_n(fn);
    bond -> setForce_st(fs);
}
