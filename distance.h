double distance(double *pos1, double *pos2)
{
  double D;
  D = sqrt(pow(pos2[0]-pos1[0],2.)+pow(pos2[1]-pos1[1],2.)+pow(pos2[2]-pos1[2],2.));
  return D;
}
