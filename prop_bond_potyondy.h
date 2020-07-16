class Bond
{
public:
  Bond() :
    sigmac_(0.),
    tauc_(0.),
    rad_(0.),
    kn_(0.),
    ks_(0.),
    partid_{0, 0},
    inertia_{0., 0.},
    force_n_{0., 0., 0.},
    force_s_{0., 0., 0.},
    torque_n_{0., 0., 0.},
    torque_s_{0., 0., 0.},
    torque_st_{0., 0., 0.}
  {}

  Bond(double sigma, double tau, double radius, double kn, double ks, double area, double kb, double kt) :
    sigmac_(sigma),
    tauc_(tau),
    rad_(radius),
    kn_(kn),
    ks_(ks)
  {}

  Bond(int i, int j, double I, double J) :
    partid_{i, j},
    inertia_{I, J}
  {}

  Bond(double fnx, double fny, double fnz, double fsx, double fsy, double fsz, double tnx, double tny, double tnz, double tstx, double tsty, double tstz, double tsrx, double tsry, double tsrz, double ttrx, double ttry, double ttrz) :
    force_n_{fnx, fny, fnz},
    force_s_{fsx, fsy, fsz},
    torque_n_{tnx, tny, tnz},
    torque_s_{tstx, tsty, tstz},
    torque_st_{tsrx, tsry, tsrz}
  {}


  ~Bond() {}

  double getSigma()
  { return sigmac_; }

  double getTau()
  { return tauc_; }

  double getRad()
  { return rad_; }

  double getStiffNorm()
  { return kn_; }

  double getStiffShear()
  { return ks_; }

  int* getIds()
  { return partid_; }

  double* getInertia()
  { return inertia_; }

  double* getForce_n()
  { return force_n_; }

  double* getForce_s()
  { return force_s_; }

  double* getTorque_n()
  { return torque_n_; }

  double* getTorque_s()
  { return torque_s_; }

  double* getTorque_st()
  { return torque_st_; }

  void setSigma(double sigma)
  { sigmac_ = sigma; }

  void setTau(double tau)
  { tauc_ = tau; }

  void setRad(double radius)
  { rad_ = radius; }

  void setStiffNorm(double kn)
  { kn_ = kn; }

  void setStiffShear(double ks)
  { ks_ = ks; }

  void setIds(int* ids)
  {
    partid_[0] = ids[0];
    partid_[1] = ids[1];
  }

  void setInertia(double* in)
  {
    inertia_[0] = in[0];
    inertia_[1] = in[1];
  }

  void setForce_n(double* fn)
  {
    force_n_[0] = fn[0];
    force_n_[1] = fn[1];
    force_n_[2] = fn[2];
  }

  void setForce_s(double* fs)
  {
    force_s_[0] = fs[0];
    force_s_[1] = fs[1];
    force_s_[2] = fs[2];
  }

  void setTorque_n(double* tn)
  {
    torque_n_[0] = tn[0];
    torque_n_[1] = tn[1];
    torque_n_[2] = tn[2];
  }

  void setTorque_s(double* ts)
  {
    torque_s_[0] = ts[0];
    torque_s_[1] = ts[1];
    torque_s_[2] = ts[2];
  }

  void setTorque_st(double* tst)
  {
    torque_st_[0] = tst[0];
    torque_st_[1] = tst[1];
    torque_st_[2] = tst[2];
  }


private:
  double sigmac_, tauc_, rad_, kn_, ks_, inertia_[2];
  int partid_[2];
  double force_n_[3], force_s_[3], torque_n_[3], torque_s_[3], torque_st_[3];

};




