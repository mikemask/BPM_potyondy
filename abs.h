double abs(double *vec)
{
    double norm = sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]);

    return norm;
}
