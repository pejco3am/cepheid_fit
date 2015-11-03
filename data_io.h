void read_data(data_structure *data, char const  *filname);
void read_periods(ephemeris *eph, char const *filname, const int cep_fit);
void read_beta(double *beta, char const *filname);
void read_fc(double *fc, char const *filname, const double period);
int read_cep_pars(cep_pars *cep, char const *filname, const int id);
void write_cep_pars(cep_pars *cep, char const *filname, const int id);
void write_periods(ephemeris *eph, char const *filname, const int id);
