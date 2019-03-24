
double myinteglin(int n, double *x, double *y, double z){
    double yofz = mylinterp(n, x, y, z);

    double integral = 0; /* Linear integration by average square between the ponts */
    int i = 0;
    while (x[i+1] < z) {
        integral += 1/2.0 * (y[i + 1] + y[i])*(x[i + 1] - x[i]);
        i++;
    }  /* We now have integral up to the point before z */
 /* We will now add the last contribution up to z with value S_z */
    integral += 1/2.0 * (yofz + y[i])*(z - x[i]);
    return integral;

}
