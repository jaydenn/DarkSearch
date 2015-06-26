
int getSamplingPars(parameterList *pL, char *filename);

int writeSamplingOutput(parameterList pL);

int writeRateOutput(parameterList pL, int detj, double *Er, double *signal, double *background, int sizeData);
