	#include <stdio.h>
	#include <stdlib.h>
	#include <math.h>
	#include <time.h>

	#define  NTAB 32
	#define  IR 2836
	#define  IM3 2147483647
	#define  XDIV (1+(IM3-1)/NTAB)
	#define  RNMX (1.0-EPS)
	#define  AM1 (1.0/IM3)
	#define  IQ 127773
	#define  IA 16807
	#define  EPS 1.2e-7

	/*******************************/
	void str_setup(char str0[], char str1[], char str2[], char str3[], char str4[]);
	int  mystrcmp( char s1[], char s2[] );
	void mystrcat( char str0[], char str[] );
	long str_len(char str[]);
	void mystrcpy( char str0[], char str[] );
	double  mystrtolf( char S[] );
	void myltos(long N, char istr[]);
	void mydtos(double N, char istr[], int decimals);
	/****************************************/
	double  ran1( long *idumPtr );
	long UniformRndInt(long jlo, long jhi, double x);
	/****************************************/
	char **cmatrix(long nrl, long nrh, long ncl, long nch);
	void free_cmatrix(char **m, long nrl, long nrh, long ncl, long nch);
	char *cvector(long nl, long nh);
	void free_cvector(char *v, long nl, long nh);
	long *lvector(long nl, long nh);
	void free_lvector(long *v, long nl, long nh);
	double *dvector(long nl, long nh);
	void free_dvector(double *v, long nl, long nh);
	double **dmatrix(long nrl, long nrh, long ncl, long nch);
	void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);
	long **lmatrix(long nrl, long nrh, long ncl, long nch);
	void free_lmatrix(long **m, long nrl, long nrh, long ncl, long nch);
	int  *ivector(long nl, long nh);
	void free_ivector(int *v, long nl, long nh);
	unsigned long *ulvector(long nl, long nh);
	void free_ulvector(unsigned long *v, long nl, long nh);
	/*************************/
	void Deaths(int model,long nloci,long *idumPtr,long age_classes,long popN,double *survival,
									long **mut_mort,long *age_of_expression,double mut_mort_effect,long *age,long *deaths);
	void Parents(int model,long *idumPtr,long popN,long nloci,long age_classes,long female,long male,double **mate,
										long *female_age, long *male_age, double *female_fertility, double *male_fertility,
										long **female_mut_fert,long **male_mut_fert,long *age_of_expression, double mut_fert_effect,
										long *motherPtr, long *fatherPtr);
	double GetSurvival(long nloci,long age_classes,long *age_of_expression,double mut_mort_effect,
																long age,double *survival,long *mut_mort);
	double GetFertility(long nloci,long age_classes,long *age_of_expression,double mut_fert_effect,
																	long age,double *fertility,long *mut_fert);
	void AgeDist(long popN, long age_classes,long *age_dist, long *age);
	void MutDist(long popN, long age_classes, long *age,long *mut_dist, long **mut);
	void FertDist(char sex,long nloci,long popN, long age_classes, long *age,double *fert_dist,
											long **mut,double *fertility,long *age_of_expression, double mut_effect);
	void MortDist(char sex,long nloci,long popN, long age_classes, long *age,double *mort_dist,
											long **mut,double *survival,long *age_of_expression, double mut_effect);
	long GetMut(long *idumPtr,long nloci,long *mut,double mu);
	void Transmission(long nloci,long *male_mut,long *female_mut,long *offspring_mut,long *idumPtr);
	void LociFixed(long nloci,long popN,long **mut,long *fixed);
	void MeanVar(long n, double *data, double *meanPtr, double *svarPtr);
	/*******/
	int main (int argc, const char * argv[]) {
					time_t t1;
					long age_classes,age,a,age_interval,m,n;
					long *age_of_expression,aoe;
					long i,tc,ti,t,nloci,nl,interval;
					long *female_age,*init_female_age,mother,*mothers,*female_fixed_mort,*female_fixed_fert;
					long *male_age,*init_male_age,father,*fathers,*male_fixed_mort,*male_fixed_fert;
					long *female_deaths,*female_age_dist,**female_age_Rdist,**female_mut_mort_Rdist,**female_mut_mort_Ldist;
					long *male_deaths,*male_age_dist,**male_age_Rdist,**male_mut_mort_Rdist,**male_mut_mort_Ldist;
					long *male_mut_mort_dist,*male_mut_fert_dist;
					long *female_mut_mort_dist,*female_mut_fert_dist;
					long **parents_age_dist;
					long idum = -175566,*cum_new_mut_mort,new_mut_mort;
					long repeats,r,popN;
					long time_count,total_time_counts;

	/* number of mortality-causing mutations carried by i-individual */
					long **female_mut_mort,**init_female_mut_mort,**female_mut_fert,**init_female_mut_fert;
					long **male_mut_mort,**init_male_mut_mort,**male_mut_fert,**init_male_mut_fert;

					double *male_fertility,*female_fertility;
					double *female_survival,**female_fert_Rdist,**female_mort_Rdist,**female_fixed_mort_dynamic,**female_mut_mort_dynamic;
					double *male_survival,**male_fert_Rdist,**male_mort_Rdist,**male_fixed_mort_dynamic,**male_mut_mort_dynamic;
					double **mate;
					double x,mu_mort,mut_mort_effect,mu_fert,mut_fert_effect,mean,var;

					char cntlfile[60],run[20],inputfile[240];
					char str[240],logpath[240],globalpath[240],agedistpath[240],fertdistpath[240],mortdistpath[240];
					char deathdistpath[240],mut_mortdistpath[240],mut_fertdistpath[240];
					char birthdistpath[240],outpath_time[240],parentspath[240],mut_mort_locipath[240];

					FILE *input,*input_cntl,*logout,*agedistout,*mut_mortdistout,*mut_fertdistout,*mut_mort_locifile;
					FILE *parentsdistout,*outfile_time,*fertdistout,*deathdistout,*birthdistout,*mortdistout;
	/*********************************************/
					mystrcpy(cntlfile,argv[1]);
					/*printf("\nRead instructions from %s\n",cntlfile);*/
					input_cntl=fopen(cntlfile,"r");
					if(input_cntl == NULL) {printf("can't open %s\n",cntlfile); exit(1); }
	/* Read global path */
					if(fscanf(input_cntl,"%s %s ",str,globalpath) != 2) {printf("can't read global path\n"); exit(1); }
					if(mystrcmp(str,"Global") != 1) {printf("can't read Global [%s]\n",str); exit(1); }
					printf("Global path:\t%s\n",globalpath);
	/* read run */
					if(fscanf(input_cntl,"%s %s ",str,run) != 2) {printf("can't read run\n"); exit(1); }
					if(mystrcmp(str,"Run") != 1) {printf("can't read Run [%s]\n",str); exit(1); }
					mystrcpy(logpath,globalpath);
					mystrcat(logpath,run);
					mystrcat(logpath,"/");
					mystrcat(logpath,run);
					mystrcat(logpath,"_log.txt");
					printf("Log path:\t%s\n",logpath);
					logout=fopen(logpath,"w");
					if(logout == NULL) {printf("Can't open %s\n",logpath); exit(1); }
					t1 = time(NULL);
					fprintf(logout,"Reproduce\n%s\n",ctime(&t1));
					fprintf(logout,"Global path\t%s\n",globalpath);
					fprintf(logout,"Run\t%s\n",run);
	/* read number of age classes */
					if(fscanf(input_cntl,"%s %ld ",str,&age_classes) != 2) {printf("can't read age_classes\n"); exit(1); }
					if(mystrcmp(str,"Age_Classes") != 1) {printf("can't read Age_Classes [%s]\n",str); exit(1); }
					printf("Age_Classes:\t%ld\n",age_classes);
					fprintf(logout,"Age Classes\t%ld\n",age_classes);
	/* read age class interval for output */
					if(fscanf(input_cntl,"%s %ld ",str,&age_interval) != 2) {printf("can't read age interval\n"); exit(1); }
					if(mystrcmp(str,"Age_class_interval") != 1) {printf("can't read Age_class_interval [%s]\n",str); exit(1); }
					printf("Age Interval\t%ld\n",age_interval);
					fprintf(logout,"Age Interval\t%ld\n",age_interval);
	/* read population (female = maless) size */
					if(fscanf(input_cntl,"%s %ld ",str,&popN) != 2) {printf("can't read popN\n"); exit(1); }
					if(mystrcmp(str,"Population_size") != 1) {printf("can't read Population_size [%s]\n",str); exit(1); }
					printf("Population Size\t%ld\n",popN);
					fprintf(logout,"Population Size\t%ld\n",popN);
	/* read mut_mort mutation rate */
					if(fscanf(input_cntl,"%s %le ",str,&mu_mort) != 2) {printf("can't read mu_mort\n"); exit(1); }
					if(mystrcmp(str,"mort_mutation_rate") != 1) {printf("can't read mort_mutation_rate [%s]\n",str); exit(1); }
					printf("Mortality-Causing Mutation Rate\t%lf\n",mu_mort);
	/* read mort mutation effect */
					if(fscanf(input_cntl,"%s %lf ",str,&mut_mort_effect) != 2) {printf("can't read mut_mort_effect\n"); exit(1); }
					if(mystrcmp(str,"mort_mutation_effect") != 1) {printf("can't read mort_mutation_effect [%s]\n",str); exit(1); }
					printf("Mortality-Causing Mutation Effect\t%lf\n",mut_mort_effect);
	/* read mut_fert mutation rate */
					if(fscanf(input_cntl,"%s %le ",str,&mu_fert) != 2) {printf("can't read mu_fert\n"); exit(1); }
					if(mystrcmp(str,"fert_mutation_rate") != 1) {printf("can't read fert_mutation_rate [%s]\n",str); exit(1); }
					printf("Fertility-Affecting Mutation Rate\t%lf\n",mu_fert);
					fprintf(logout,"Mortality-Causing Mutation Rate\t%lf\n",mu_mort);
					fprintf(logout,"Fertility-Affecting Mutation Rate\t%lf\n",mu_fert);
	/* read fert mutation effect */
					if(fscanf(input_cntl,"%s %lf ",str,&mut_fert_effect) != 2) {printf("can't read mut_fert_effect\n"); exit(1); }
					if(mystrcmp(str,"fert_mutation_effect") != 1) {printf("can't read fert_mutation_effect [%s]\n",str); exit(1); }
					printf("Fertility-Affecting Mutation Effect\t%lf\n",mut_fert_effect);
					fprintf(logout,"Mortality-Causing Mutation Effect\t%lf\n",mut_mort_effect);
					fprintf(logout,"Fertility-Affecting Mutation Effect\t%lf\n",mut_fert_effect);
	/* read age of expression */
					if(fscanf(input_cntl,"%s %ld ",str,&aoe) != 2) {printf("can't read aoe\n"); exit(1); }
					if(mystrcmp(str,"Age_of_expression") != 1) {printf("can't read Age_of_expression [%s]\n",str); exit(1); }
					printf("Initial Age of Expression\t%ld\n",aoe);
					fprintf(logout,"Initial Age of Expression\t%ld\n",aoe);
	/* read number of loci */
					if(fscanf(input_cntl,"%s %ld ",str,&nloci) != 2) {printf("can't read number of loci\n"); exit(1); }
					if(mystrcmp(str,"Number_of_loci") != 1) {printf("can't read Number_of_loci [%s]\n",str); exit(1); }
					printf("Number of loci\t%ld\n",nloci);
					fprintf(logout,"Number of loci\t%ld\n",nloci);
					fprintf(logout,"\n");
					age_of_expression=lvector(0,nloci);
					age_of_expression[0]=0;
					if(nloci >= 1) {
													age_of_expression[1]=aoe;
													if(nloci > 1) for(a=2; a<=nloci; a++) {
																													age_of_expression[a]=age_of_expression[a-1]+1;
																					}
					}
					fprintf(logout,"locus\tage_of_expression\n");
					if(nloci <= 0) fprintf(logout,"no mutations\n");
					else { for(i=1; i<=nloci; i++) {
																				fprintf(logout,"%ld\t%ld\n",i,age_of_expression[i]);
												} }
					if(nloci <= 5) interval=1;
					if(nloci > 5 && nloci <= 10) interval=2;
					if(nloci > 10) interval=5;
	/* read male survival file */
					male_survival=dvector(0,age_classes);
					if(fscanf(input_cntl,"%s %s ",str,inputfile) != 2) {printf("can't read male survival file\n"); exit(1); }
					if(mystrcmp(str,"Male_Survival") != 1) {printf("can't read Male_Survival [%s]\n",str); exit(1); }
					printf("Male Survival file\t%s\n",inputfile);
					input=fopen(inputfile,"r");
					if(input == NULL) {printf("Can't open %s\n",inputfile); exit(1); }
					if(fscanf(input,"%s %ld ",str,&a) != 2) {printf("can't read header initial distribution file\n"); exit(1); }
					if(mystrcmp(str,"Age_Classes") != 1) {printf("can't read Age_Classes [%s]\n",str); exit(1); }
					if(a != age_classes) {printf("Age classes [%ld] != %ld\n",age_classes,a); exit(1); }
					for(age=0; age<=age_classes; age++) {
													if(fscanf(input,"%ld %lf ",&a,&male_survival[age]) != 2) {printf("can't read male survival age = %ld\n",age); exit(1); }
													if(a != age) {printf("Age_class [%ld] != %ld\n",age,a); exit(1); }
													if(male_survival[age] < 0.0 || male_survival[age] > 1.0) {printf("male survival[%ld] = %lf out of range\n",age,male_survival[age]); exit(1); }
													/*printf("age=%ld\tsurvival = %lf\n",age,survival[age]);*/
					}
					fclose(input);
	/* read female survival file */
					female_survival=dvector(0,age_classes);
					if(fscanf(input_cntl,"%s %s ",str,inputfile) != 2) {printf("can't read female survival file\n"); exit(1); }
					if(mystrcmp(str,"Female_Survival") != 1) {printf("can't read Female_Survival [%s]\n",str); exit(1); }
					printf("Female Survival file\t%s\n",inputfile);
					input=fopen(inputfile,"r");
					if(input == NULL) {printf("Can't open %s\n",inputfile); exit(1); }
					if(fscanf(input,"%s %ld ",str,&a) != 2) {printf("can't read header initial distribution file\n"); exit(1); }
					if(mystrcmp(str,"Age_Classes") != 1) {printf("can't read Age_Classes [%s]\n",str); exit(1); }
					if(a != age_classes) {printf("Age classes [%ld] != %ld\n",age_classes,a); exit(1); }
					for(age=0; age<=age_classes; age++) {
													if(fscanf(input,"%ld %lf ",&a,&female_survival[age]) != 2) {printf("can't read female survival age = %ld\n",age); exit(1); }
													if(a != age) {printf("Age_class [%ld] != %ld\n",age,a); exit(1); }
													if(female_survival[age] < 0.0 || female_survival[age] > 1.0) {printf("female survival[%ld] = %lf out of range\n",age,female_survival[age]); exit(1); }
													/*printf("age=%ld\tsurvival = %lf\n",age,survival[age]);*/
					}
					fclose(input);
	/* read mating file */
					mate=dmatrix(0,age_classes,0,age_classes);
					if(fscanf(input_cntl,"%s %s ",str,inputfile) != 2) {printf("can't read mating file\n"); exit(1); }
					if(mystrcmp(str,"Mate") != 1) {printf("can't read Mate [%s]\n",str); exit(1); }
					printf("Mating file\t%s\n",inputfile);
					input=fopen(inputfile,"r");
					if(input == NULL) {printf("Can't open %s\n",inputfile); exit(1); }
					if(fscanf(input,"%s %ld ",str,&a) != 2) {printf("can't read header mating file\n"); exit(1); }
					if(mystrcmp(str,"Age_Classes") != 1) {printf("can't read Age_Classes [%s]\n",str); exit(1); }
					if(a != age_classes) {printf("Age classes [%ld] != %ld\n",age_classes,a); exit(1); }
					for(age=0; age<=age_classes; age++) { /* age is male age */
													if(fscanf(input,"%ld ",&a) != 1) {printf("can't read male age class = %ld\n",age); exit(1); }
													if(a != age) {printf("Age_class [%ld] != %ld\n",age,a); exit(1); }
													for(a=0; a<=age_classes; a++) { /* a is female age */
																					if(fscanf(input,"%lf ",&mate[age][a]) != 1) {printf("can't read mate age = %ld\n",a); exit(1); }
																					if(mate[age][a] < 0.0 || mate[age][a] > 1.0) {printf("mate[%ld][%ld] = %lf out of range\n",age,a,mate[age][a]); exit(1); }
													}
					}
					fclose(input);
	/* read male birth rate file */
					male_fertility=dvector(0,age_classes);
					if(fscanf(input_cntl,"%s %s ",str,inputfile) != 2) {printf("can't read male birth rate file\n"); exit(1); }
					if(mystrcmp(str,"Male_fertility") != 1) {printf("can't read Male_fertility [%s]\n",str); exit(1); }
					printf("Male Fertility Rate file\t%s\n",inputfile);
					input=fopen(inputfile,"r");
					if(input == NULL) {printf("Can't open %s\n",inputfile); exit(1); }
					if(fscanf(input,"%s %ld ",str,&a) != 2) {printf("can't read header male birth rate file\n"); exit(1); }
					if(mystrcmp(str,"Age_Classes") != 1) {printf("can't read Age_Classes [%s]\n",str); exit(1); }
					if(a != age_classes) {printf("Age classes [%ld] != %ld\n",age_classes,a); exit(1); }
					for(age=0; age<=age_classes; age++) {
													if(fscanf(input,"%ld %lf ",&a,&male_fertility[age]) != 2) {printf("can't read male birth age = %ld\n",age); exit(1); }
													if(a != age) {printf("Age_class [%ld] != %ld\n",age,a); exit(1); }
													if(male_fertility[age] < 0.0 || male_fertility[age] > 1.0) {printf("male_fertility[%ld] = %lf out of range\n",age,male_fertility[age]); exit(1); }
													/* printf("age=%ld\tbirth = %lf\n",age,birth[age]);*/
					}
					fclose(input);
	/* read female birth rate file */
					female_fertility=dvector(0,age_classes);
					if(fscanf(input_cntl,"%s %s ",str,inputfile) != 2) {printf("can't read female birth rate file\n"); exit(1); }
					if(mystrcmp(str,"Female_fertility") != 1) {printf("can't read Female_fertility [%s]\n",str); exit(1); }
					printf("Female Fertility Rate file\t%s\n",inputfile);
					input=fopen(inputfile,"r");
					if(input == NULL) {printf("Can't open %s\n",inputfile); exit(1); }
					if(fscanf(input,"%s %ld ",str,&a) != 2) {printf("can't read header female birth rate file\n"); exit(1); }
					if(mystrcmp(str,"Age_Classes") != 1) {printf("can't read Age_Classes [%s]\n",str); exit(1); }
					if(a != age_classes) {printf("Age classes [%ld] != %ld\n",age_classes,a); exit(1); }
					for(age=0; age<=age_classes; age++) {
													if(fscanf(input,"%ld %lf ",&a,&female_fertility[age]) != 2) {printf("can't read female birth age = %ld\n",age); exit(1); }
													if(a != age) {printf("Age_class [%ld] != %ld\n",age,a); exit(1); }
													if(female_fertility[age] < 0.0 || female_fertility[age] > 1.0) {printf("female_fertility[%ld] = %lf out of range\n",age,female_fertility[age]); exit(1); }
													/* printf("age=%ld\tbirth = %lf\n",age,birth[age]);*/
					}
					fclose(input);
	/* read repeats of simulation */
					if(fscanf(input_cntl,"%s %ld ",str,&repeats) != 2) {printf("can't read repeats\n"); exit(1); }
					if(mystrcmp(str,"Repeats") != 1) {printf("can't read Repeats [%s]\n",str); exit(1); }
					printf("repeats = %ld\n",repeats);
	/* read time_count of simulation */
					if(fscanf(input_cntl,"%s %ld ",str,&time_count) != 2) {printf("can't read time_count\n"); exit(1); }
					if(mystrcmp(str,"Time_count") != 1) {printf("can't read Time_count [%s]\n",str); exit(1); }
					printf("time_count = %ld\n",time_count);
	/* read total_time_counts of simulation */
					if(fscanf(input_cntl,"%s %ld ",str,&total_time_counts) != 2) {printf("can't read total_time_counts\n"); exit(1); }
					if(mystrcmp(str,"Total_time_count") != 1) {printf("can't read Total_time_count [%s]\n",str); exit(1); }
					printf("Total_time_count = %ld\n",total_time_counts);
					fclose(input_cntl);

	/** arrays */
					init_male_mut_mort=lmatrix(1,popN,0,nloci);
					male_mut_mort=lmatrix(1,popN,0,nloci); /*male mut_mortality by locus [0] holds total */
					init_male_mut_fert=lmatrix(1,popN,0,nloci);
					male_mut_fert=lmatrix(1,popN,0,nloci); /*male mut_fertility by locus [0] holds total */
					male_mut_mort_dist=lvector(0,age_classes);
					male_mut_fert_dist=lvector(0,age_classes);
					male_deaths=lvector(0,age_classes);
					fathers=lvector(0,age_classes);
					init_male_age=lvector(1,popN);
					male_age=lvector(1,popN); /*male age */
					male_age_dist=lvector(0,age_classes);

					male_age_Rdist=lmatrix(0,repeats,0,age_classes);
					male_fert_Rdist=dmatrix(0,repeats,0,age_classes);
					male_mort_Rdist=dmatrix(0,repeats,0,age_classes);
					male_mut_mort_Rdist=lmatrix(0,repeats,0,age_classes);
					male_mut_mort_Ldist=lmatrix(0,repeats,0,nloci); /* number of mutations at each locus in pop; [0]=total */

					male_fixed_mort=lvector(0,nloci);
					male_fixed_fert=lvector(0,nloci);
					male_fixed_mort_dynamic=dmatrix(0,total_time_counts,0,repeats);
					male_mut_mort_dynamic=dmatrix(0,total_time_counts,0,repeats);

					init_female_mut_mort=lmatrix(1,popN,0,nloci);
					female_mut_mort=lmatrix(1,popN,0,nloci); /*female mut_mortality by locus [0] holds total */
					init_female_mut_fert=lmatrix(1,popN,0,nloci);
					female_mut_fert=lmatrix(1,popN,0,nloci); /*female mut_fertility by locus [0] holds total */
					female_mut_mort_dist=lvector(0,age_classes);
					female_mut_fert_dist=lvector(0,age_classes);
					female_deaths=lvector(0,age_classes);
					mothers=lvector(0,age_classes);
					female_age_dist=lvector(0,age_classes);

					female_fert_Rdist=dmatrix(0,repeats,0,age_classes);
					female_mort_Rdist=dmatrix(0,repeats,0,age_classes);
					female_mut_mort_Rdist=lmatrix(0,repeats,0,age_classes);
					female_age_Rdist=lmatrix(0,repeats,0,age_classes);
					female_mut_mort_Ldist=lmatrix(0,repeats,0,nloci); /* number of mutations at each locus in pop; [0]=total */

					init_female_age=lvector(1,popN);
					female_age=lvector(1,popN); /*female age */
					female_fixed_mort=lvector(0,nloci);
					female_fixed_fert=lvector(0,nloci);
					female_fixed_mort_dynamic=dmatrix(0,total_time_counts,0,repeats);
					female_mut_mort_dynamic=dmatrix(0,total_time_counts,0,repeats);

					cum_new_mut_mort=lvector(1,repeats);
					parents_age_dist=lmatrix(0,age_classes,0,age_classes);

	/* male survival*/
					fprintf(logout,"\nMale Survival\tMutations\n");
					fprintf(logout,"Age Class\t");
					for(m=1; m<=nloci; m++) for(i=0; i<=2; i++) {fprintf(logout,"Locus #%ld mut=%ld\t",m,i); }
					fprintf(logout,"\n");
					for(age=1; age<=age_classes; age++) {
													fprintf(logout,"%ld\t",age);
													for(m=1; m<=nloci; m++) {
																					for(n=0; n<=nloci; n++) {male_mut_mort[1][n]=0; }
																					for(i=0; i<=2; i++) {
																													if(i == 0) {male_mut_mort[1][m]=0; male_mut_mort[1][0]=0; }
																													if(i == 1) {male_mut_mort[1][m]=1; male_mut_mort[1][0]=1; }
																													if(i == 2) {male_mut_mort[1][m]=2; male_mut_mort[1][0]=2; }
																													x=GetSurvival(nloci,age_classes,age_of_expression,mut_mort_effect,age,male_survival,male_mut_mort[1]);
																													fprintf(logout,"%lf\t",x);
																					}
													}
													fprintf(logout,"\n");
													/*{printf("pause..");getchar();printf("\n");}*/
					}
	/* female survival*/
					fprintf(logout,"\nFemale Survival\tMutations\n");
					fprintf(logout,"Age Class\t");
					for(m=1; m<=nloci; m++) for(i=0; i<=2; i++) {fprintf(logout,"Locus #%ld mut=%ld\t",m,i); }
					fprintf(logout,"\n");
					for(age=1; age<=age_classes; age++) {
													fprintf(logout,"%ld\t",age);
													for(m=1; m<=nloci; m++) {
																					for(n=0; n<=nloci; n++) {female_mut_mort[1][n]=0; }
																					for(i=0; i<=2; i++) {
																													if(i == 0) {female_mut_mort[1][m]=0; female_mut_mort[1][0]=0; }
																													if(i == 1) {female_mut_mort[1][m]=1; female_mut_mort[1][0]=1; }
																													if(i == 2) {female_mut_mort[1][m]=2; female_mut_mort[1][0]=2; }
																													x=GetSurvival(nloci,age_classes,age_of_expression,mut_mort_effect,age,female_survival,female_mut_mort[1]);
																													fprintf(logout,"%lf\t",x);
																					}
													}
													fprintf(logout,"\n");
					}
	/* male fertility*/
					fprintf(logout,"\nMale Fertility\tMutations\n");
					fprintf(logout,"Age Class\t");
					for(m=1; m<=nloci; m++) for(i=0; i<=2; i++) {fprintf(logout,"Locus #%ld mut=%ld\t",m,i); }
					fprintf(logout,"\n");
					for(age=1; age<=age_classes; age++) {
													fprintf(logout,"%ld\t",age);
													for(m=1; m<=nloci; m++) {
																					for(n=0; n<=nloci; n++) {male_mut_fert[1][n]=0; }
																					for(i=0; i<=2; i++) {
																													if(i == 0) {male_mut_fert[1][m]=0; male_mut_fert[1][0]=0; }
																													if(i == 1) {male_mut_fert[1][m]=1; male_mut_fert[1][0]=1; }
																													if(i == 2) {male_mut_fert[1][m]=2; male_mut_fert[1][0]=2; }
																													x=GetFertility(nloci,age_classes,age_of_expression,mut_fert_effect,age,male_fertility,male_mut_fert[1]);
																													fprintf(logout,"%lf\t",x);
																					}
													}
													fprintf(logout,"\n");
													/*{printf("pause..");getchar();printf("\n");}*/
					}
	/* female fertility*/
					fprintf(logout,"\nFemale Fertility\tMutations\n");
					fprintf(logout,"Age Class\t");
					for(m=1; m<=nloci; m++) for(i=0; i<=2; i++) {fprintf(logout,"Locus #%ld mut=%ld\t",m,i); }
					fprintf(logout,"\n");
					for(age=1; age<=age_classes; age++) {
													fprintf(logout,"%ld\t",age);
													for(m=1; m<=nloci; m++) {
																					for(n=0; n<=nloci; n++) {female_mut_fert[1][n]=0; }
																					for(i=0; i<=2; i++) {
																													if(i == 0) {female_mut_fert[1][m]=0; female_mut_fert[1][0]=0; }
																													if(i == 1) {female_mut_fert[1][m]=1; female_mut_fert[1][0]=1; }
																													if(i == 2) {female_mut_fert[1][m]=2; female_mut_fert[1][0]=2; }
																													x=GetFertility(nloci,age_classes,age_of_expression,mut_fert_effect,age,female_fertility,female_mut_fert[1]);
																													fprintf(logout,"%lf\t",x);
																					}
													}
													fprintf(logout,"\n");
													/*{printf("pause..");getchar();printf("\n");}*/
					}

					fprintf(logout,"\nMate Matrix\n\t0\t");
					for(a=1; a<=age_classes; a++) {fprintf(logout,"%ld to %ld-\t",(a-1)*age_interval,(a)*age_interval); }
					fprintf(logout,"\n");
					fprintf(logout,"0\t");
					for(a=1; a<=age_classes; a++) {fprintf(logout,"%lf\t",mate[0][a]); }
					fprintf(logout,"\n");
					for(age=1; age<=age_classes; age++) {
													fprintf(logout,"%ld to %ld-\t",(age-1)*age_interval,(age)*age_interval);
													for(a=1; a<=age_classes; a++) {fprintf(logout,"%lf\t",mate[age][a]); }
													fprintf(logout,"\n");
					}
					fprintf(logout,"End input data\n\n");
					fclose(logout);

	/* file for age distributions */
					mystrcpy(agedistpath,globalpath);
					mystrcat(agedistpath,run);
					mystrcat(agedistpath,"/");
					mystrcat(agedistpath,run);
					mystrcat(agedistpath,"_agedist.txt");
					agedistout=fopen(agedistpath,"w");
					if(agedistout == NULL) {printf("Can't open %s\n",agedistpath); exit(1); }
					fprintf(agedistout,"Age-Specific distribution of population size\n");
					fprintf(agedistout,"popN\t%ld\nnloci\t%ld\nmut_mort_rate\t%lf\nmut_mort_effect\t%lf\nmut_fert_rate\t%lf\nmut_fert_effect\t%lf\n",popN,nloci,mu_mort,mut_mort_effect,mu_fert,mut_fert_effect);
					fprintf(agedistout,"Repeat (at end)\t");
					fprintf(agedistout,"Total (male)\t");
					for(age=1; age<=(age_classes-1); age++) fprintf(agedistout,"(male) %ld to %ld-\t",(age-1)*age_interval,(age)*age_interval);
					fprintf(agedistout,"(male) >%ld\t",age_classes*age_interval);
					fprintf(agedistout,"Total (female)\t");
					for(age=1; age<=(age_classes-1); age++) fprintf(agedistout,"(female) %ld to %ld-\t",(age-1)*age_interval,(age)*age_interval);
					fprintf(agedistout,"(female) >%ld\n",age_classes*age_interval);
					fclose(agedistout);

	/* file for mut_mortality distributions */
					mystrcpy(mut_mortdistpath,globalpath);
					mystrcat(mut_mortdistpath,run);
					mystrcat(mut_mortdistpath,"/");
					mystrcat(mut_mortdistpath,run);
					mystrcat(mut_mortdistpath,"_mut_mortdist.txt");
					mut_mortdistout=fopen(mut_mortdistpath,"w");
					if(mut_mortdistout == NULL) {printf("Can't open %s\n",mut_mortdistpath); exit(1); }
					fprintf(mut_mortdistout,"Age-Specific distribution of mortality-affecting\n");
					fprintf(mut_mortdistout,"popN\t%ld\nnloci\t%ld\nmut_mort_rate\t%lf\nmut_mort_effect\t%lf\nmut_fert_rate\t%lf\nmut_fert_effect\t%lf\n",popN,nloci,mu_mort,mut_mort_effect,mu_fert,mut_fert_effect);
					fprintf(mut_mortdistout,"Repeat (at end)\t");
					fprintf(mut_mortdistout,"Loci fixed male out of %ld\tTotal (male)\t",nloci*popN);
					for(age=1; age<=(age_classes-1); age++) fprintf(mut_mortdistout,"(male) %ld to %ld-\t",(age-1)*age_interval,(age)*age_interval);
					fprintf(mut_mortdistout,"(male) >%ld\t",age_classes*age_interval);
					fprintf(mut_mortdistout,"Loci fixed female out of %ld\tTotal (female)\t",nloci*popN);
					for(age=1; age<=(age_classes-1); age++) fprintf(mut_mortdistout,"(female) %ld to %ld-\t",(age-1)*age_interval,(age)*age_interval);
					fprintf(mut_mortdistout,"(female) >%ld\n",age_classes*age_interval);
					fclose(mut_mortdistout);

	/* file for mut_fertility distributions */
					mystrcpy(mut_fertdistpath,globalpath);
					mystrcat(mut_fertdistpath,run);
					mystrcat(mut_fertdistpath,"/");
					mystrcat(mut_fertdistpath,run);
					mystrcat(mut_fertdistpath,"_mut_fertdist.txt");
					mut_fertdistout=fopen(mut_fertdistpath,"w");
					if(mut_fertdistout == NULL) {printf("Can't open %s\n",mut_fertdistpath); exit(1); }
					fprintf(mut_fertdistout,"Age-Specific distribution of fertility-affecting\n");
					fprintf(mut_fertdistout,"popN\t%ld\nnloci\t%ld\nmut_mort_rate\t%lf\nmut_mort_effect\t%lf\nmut_fert_rate\t%lf\nmut_fert_effect\t%lf\n",popN,nloci,mu_mort,mut_mort_effect,mu_fert,mut_fert_effect);
					fprintf(mut_fertdistout,"Repeat (at end)\t");
					fprintf(mut_fertdistout,"Loci fixed male out of %ld\tTotal (male)\t",nloci*popN);
					for(age=1; age<=(age_classes-1); age++) fprintf(mut_fertdistout,"(male) %ld to %ld-\t",(age-1)*age_interval,(age)*age_interval);
					fprintf(mut_fertdistout,"(male) >%ld\t",age_classes*age_interval);
					fprintf(mut_fertdistout,"Loci fixed female out of %ld\tTotal (female)\t",nloci*popN);
					for(age=1; age<=(age_classes-1); age++) fprintf(mut_fertdistout,"(female) %ld to %ld-\t",(age-1)*age_interval,(age)*age_interval);
					fprintf(mut_fertdistout,"(female) >%ld\n",age_classes*age_interval);
					fclose(mut_fertdistout);

	/* file for death distributions */
					mystrcpy(deathdistpath,globalpath);
					mystrcat(deathdistpath,run);
					mystrcat(deathdistpath,"/");
					mystrcat(deathdistpath,run);
					mystrcat(deathdistpath,"_deathsdist.txt");
					deathdistout=fopen(deathdistpath,"w");
					if(deathdistout == NULL) {printf("Can't open %s\n",deathdistpath); exit(1); }
					fprintf(deathdistout,"Age-Specific distribution of deaths\n");
					fprintf(deathdistout,"popN\t%ld\nnloci\t%ld\nmut_mort_rate\t%lf\nmut_mort_effect\t%lf\nmut_fert_rate\t%lf\nmut_fert_effect\t%lf\n",popN,nloci,mu_mort,mut_mort_effect,mu_fert,mut_fert_effect);
					fprintf(deathdistout,"Repeat (at end)\t");
					fprintf(deathdistout,"Total (male)\t");
					for(age=1; age<=(age_classes-1); age++) fprintf(deathdistout,"(male) %ld to %ld-\t",(age-1)*age_interval,(age)*age_interval);
					fprintf(deathdistout,"(male) >%ld\t",age_classes*age_interval);
					fprintf(deathdistout,"Total (female)\t");
					for(age=1; age<=(age_classes-1); age++) fprintf(deathdistout,"(female) %ld to %ld-\t",(age-1)*age_interval,(age)*age_interval);
					fprintf(deathdistout,"(female) >%ld\n",age_classes*age_interval);
					fclose(deathdistout);

	/* file for birth rate distributions */
					mystrcpy(birthdistpath,globalpath);
					mystrcat(birthdistpath,run);
					mystrcat(birthdistpath,"/");
					mystrcat(birthdistpath,run);
					mystrcat(birthdistpath,"_birthsdist.txt");
					birthdistout=fopen(birthdistpath,"w");
					if(birthdistout == NULL) {printf("Can't open %s\n",birthdistpath); exit(1); }
					fprintf(birthdistout,"Age-Specific distribution of births by father's and mother's age\n");
					fprintf(birthdistout,"popN\t%ld\nnloci\t%ld\nmut_mort_rate\t%lf\nmut_mort_effect\t%lf\nmut_fert_rate\t%lf\nmut_fert_effect\t%lf\n",popN,nloci,mu_mort,mut_mort_effect,mu_fert,mut_fert_effect);
					fprintf(birthdistout,"Repeat (at end)\t");
					fprintf(birthdistout,"Total (male)\t");
					for(age=1; age<=(age_classes-1); age++) fprintf(birthdistout,"(male) %ld to %ld-\t",(age-1)*age_interval,(age)*age_interval);
					fprintf(birthdistout,"(male) >%ld\t",age_classes*age_interval);
					fprintf(birthdistout,"Total (female)\t");
					for(age=1; age<=(age_classes-1); age++) fprintf(birthdistout,"(female) %ld to %ld-\t",(age-1)*age_interval,(age)*age_interval);
					fprintf(birthdistout,"(female) >%ld\n",age_classes*age_interval);
					fclose(birthdistout);

	/* file for fertility rate distributions */
					mystrcpy(fertdistpath,globalpath);
					mystrcat(fertdistpath,run);
					mystrcat(fertdistpath,"/");
					mystrcat(fertdistpath,run);
					mystrcat(fertdistpath,"_fertdist.txt");
					fertdistout=fopen(fertdistpath,"w");
					if(fertdistout == NULL) {printf("Can't open %s\n",fertdistpath); exit(1); }
					fprintf(fertdistout,"Average Age-Specific fertility by age\n");
					fprintf(fertdistout,"popN\t%ld\nnloci\t%ld\nmut_mort_rate\t%lf\nmut_mort_effect\t%lf\nmut_fert_rate\t%lf\nmut_fert_effect\t%lf\n",popN,nloci,mu_mort,mut_mort_effect,mu_fert,mut_fert_effect);
					fprintf(fertdistout,"Repeat (average fertility over individuals taken at end)\t");
					fprintf(fertdistout,"Fertility (male, total fertility in age group divided by number of individuals in age group)\t");
					for(age=1; age<=(age_classes-1); age++) fprintf(fertdistout,"(male) %ld to %ld-\t",(age-1)*age_interval,(age)*age_interval);
					fprintf(fertdistout,"(male) >%ld\t",age_classes*age_interval);
					fprintf(fertdistout,"Fertility (female, total fertility in age group divided by number of individuals in age group)\t");
					for(age=1; age<=(age_classes-1); age++) fprintf(fertdistout,"(female) %ld to %ld-\t",(age-1)*age_interval,(age)*age_interval);
					fprintf(fertdistout,"(female) >%ld\n",age_classes*age_interval);
					fclose(fertdistout);

	/* file for mortality rate distributions */
					mystrcpy(mortdistpath,globalpath);
					mystrcat(mortdistpath,run);
					mystrcat(mortdistpath,"/");
					mystrcat(mortdistpath,run);
					mystrcat(mortdistpath,"_mortdist.txt");
					mortdistout=fopen(mortdistpath,"w");
					if(mortdistout == NULL) {printf("Can't open %s\n",mortdistpath); exit(1); }
					fprintf(mortdistout,"Average Age-Specific mortality by age\n");
					fprintf(mortdistout,"popN\t%ld\nnloci\t%ld\nmut_mort_rate\t%lf\nmut_mort_effect\t%lf\nmut_fert_rate\t%lf\nmut_fert_effect\t%lf\n",popN,nloci,mu_mort,mut_mort_effect,mu_fert,mut_fert_effect);
					fprintf(mortdistout,"Repeat (average mortality over individuals taken at end)\t");
					fprintf(mortdistout,"Mortality (total male mortality in age group divided by number of male individuals in age group)\t");
					for(age=1; age<=(age_classes-1); age++) fprintf(mortdistout,"(male) %ld to %ld-\t",(age-1)*age_interval,(age)*age_interval);
					fprintf(mortdistout,"(male) >%ld\t",age_classes*age_interval);
					fprintf(mortdistout,"Mortality (total female mortality in age group divided by number of female individuals in age group)\t");
					for(age=1; age<=(age_classes-1); age++) fprintf(mortdistout,"(female) %ld to %ld-\t",(age-1)*age_interval,(age)*age_interval);
					fprintf(mortdistout,"(female) >%ld\n",age_classes*age_interval);
					fclose(mortdistout);

	/* file for parents age distribution matrix */
					mystrcpy(parentspath,globalpath);
					mystrcat(parentspath,run);
					mystrcat(parentspath,"/");
					mystrcat(parentspath,run);
					mystrcat(parentspath,"_parentsdist.txt");
					parentsdistout=fopen(parentspath,"w");
					if(parentsdistout == NULL) {printf("Can't open %s\n",parentspath); exit(1); }
					fprintf(parentsdistout,"Age-Specific distribution of parents' (mother x father) ages\n");
					fprintf(parentsdistout,"popN\t%ld\nnloci\t%ld\nmut_mort_rate\t%lf\nmut_mort_effect\t%lf\nmut_fert_rate\t%lf\nmut_fert_effect\t%lf\n",popN,nloci,mu_mort,mut_mort_effect,mu_fert,mut_fert_effect);
					fclose(parentsdistout);

	/*******************************/
					x=ran1(&idum);
					exit(1);
	/* each male/female in the population starts with random age 1 to age_classes */
	/* each male/female initially has no mut_mortations */
	/** run for 1000 time intervals to reach steady state age distribution without mut_mortation */
					for(i=1; i<=popN; i++) {
													init_female_age[i]=UniformRndInt(1, (long) age_classes,ran1(&idum));
													init_male_age[i]=UniformRndInt(1, (long) age_classes,ran1(&idum));
													for(nl=0; nl<=nloci; nl++) {init_female_mut_mort[i][nl]=0; init_male_mut_mort[i][nl]=0; }
													for(nl=0; nl<=nloci; nl++) {init_female_mut_fert[i][nl]=0; init_male_mut_fert[i][nl]=0; }
					}

	/***************** pre run burn in to make a steady state age distribution with no mut_mortations ********************************/
					logout=fopen(logpath,"a");
					if(logout == NULL) {printf("Can't open %s\n",logpath); exit(1); }
					fprintf(logout,"\nPre-run burn-in for %ld time units\n",1000);
					printf("\nPre-run burn-in for %ld time units\n",1000);
					fclose(logout);

					m=0;
					for(t=1; t<=1000; t++) {
													m++;
													if(m == 11) m=1;

	/*intialize */
													for(age=0; age<=age_classes; age++) {
																					mothers[age]=0;
																					fathers[age]=0;
																					for(a=0; a<=age_classes; a++) parents_age_dist[age][a]=0;
													}
	/* initial age distributions */
													AgeDist(popN,(long) age_classes,female_age_dist,init_female_age);
													AgeDist(popN,(long) age_classes,male_age_dist,init_male_age);

	/* determine the number of female deaths */
													Deaths(3,nloci,&idum,age_classes,popN,female_survival,init_female_mut_mort,age_of_expression,mut_mort_effect,init_female_age,female_deaths);
													/*{printf("female deaths = %ld: pause..",female_deaths[0]);getchar();printf("\n");}*/
	/* determine the number of male deaths */
													Deaths(3,nloci,&idum,age_classes,popN,male_survival,init_male_mut_mort,age_of_expression,mut_mort_effect,init_male_age,male_deaths);
													if(m == 10) printf("%ld\t%ld\t%ld\n",t,male_deaths[0],female_deaths[0]);
													/*{printf("male deaths = %ld: pause..",male_deaths[0]);getchar();printf("\n");}*/
	/* deaths are indicated by age[i] = 0 */
													/*printf("determine parents\n");*/
													for(i=1; i<=popN; i++) {
																					if(init_female_age[i] == 0) {
	/* model 0 => mut_mortations affect males and females */
																													Parents(2,&idum,popN,nloci,age_classes,i,0,mate,init_female_age,init_male_age,female_fertility,male_fertility,
																																					init_female_mut_fert,init_male_mut_fert,age_of_expression,mut_fert_effect,
																																					&mother,&father);
																													mothers[init_female_age[mother]]++;
																													mothers[0]++;
																													fathers[init_male_age[father]]++;
																													fathers[0]++;
																													parents_age_dist[init_male_age[father]][init_female_age[mother]]++;
																					}
																					if(init_male_age[i] == 0) {
																													Parents(2,&idum,popN,nloci,age_classes,0,i,mate,init_female_age,init_male_age,female_fertility,male_fertility,
																																					init_female_mut_fert,init_male_mut_fert,age_of_expression,mut_fert_effect,
																																					&mother,&father);
																													mothers[init_female_age[mother]]++;
																													mothers[0]++;
																													fathers[init_male_age[father]]++;
																													fathers[0]++;
																													parents_age_dist[init_male_age[father]][init_female_age[mother]]++;
																					}
	/* transmission of mut_mortations of parents to offspring (not here, no mut_mortations) */
	/* new mut_mortations (not here, no mut_mortations) */

													}/*for(i=1;i<=popN;i++) */
	/* update population, deaths have alread been considered, survivors move to next age group */
													for(i=1; i<=popN; i++) {
																					if(init_female_age[i] >= age_classes) {printf("female #%ld age=%ld error1\n",i,init_female_age[i]); exit(1); }
																					init_female_age[i]=init_female_age[i]+1;
																					if(init_male_age[i] >= age_classes) {printf("male #%ld age=%ld error\n",i,init_male_age[i]); exit(1); }
																					init_male_age[i]=init_male_age[i]+1;
																					/*	for(nl=0;nl<=nloci;nl++) new_female_mut_mort[i][nl]=female_mut_mort[i][nl]; */
																					/*  individual i alread has mut_mortations; if i is born in interval, it moves from age 0 to 1, but the offspring mut_mortations stille belong to i */
																					if(init_female_age[i] <= 0) {printf("female #%ld age=%ld error2\n",i,init_female_age[i]); exit(1); }
																					if(init_male_age[i] <= 0) {printf("male #%ld age=%ld error2\n",i,init_male_age[i]); exit(1); }
	/*		if(init_female_age[i] == 1) {
	printf("ofspring female #%ld has %ld mut_mortations\n",i,init_female_mut_mort[i][0]);
	}
	if(init_male_age[i] == 1) {
	printf("ofspring male #%ld has %ld mut_mortations\n",i,init_male_mut_mort[i][0]);
	}
	*/

													}/*for(i=1;i<=popN;i++)*/
	/* store initial age distributions , age-specific death rates and birth rates  first period */
													if(t == 1) {
																					logout=fopen(logpath,"a");
																					if(logout == NULL) {printf("Can't open %s\n",logpath); exit(1); }
																					fprintf(logout,"Time\t%ld\n",t);
																					fprintf(logout,"Age Class\tMales\tMale Deaths\tFathers\tFemales\tFemale Deaths\tMothers\n");
																					for(age=0; age<=age_classes; age++) {
																													fprintf(logout,"%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\n",age,male_age_dist[age],male_deaths[age],fathers[age],
																																					female_age_dist[age],female_deaths[age],mothers[age]);
																					}
																					fclose(logout);
													}
													/*{printf("pause.. t=%ld",t);getchar();printf("\n");}*/
					}/*for(t=1;t<=1000;t++) */
	/* store final age distributions , age-specific death rates and birth rates  final period */
					AgeDist(popN,(long) age_classes,female_age_dist,init_female_age);
					AgeDist(popN,(long) age_classes,male_age_dist,init_male_age);
					logout=fopen(logpath,"a");
					if(logout == NULL) {printf("Can't open %s\n",logpath); exit(1); }
					fprintf(logout,"Burn-out complete after\t%ld\nLast Interval\n",t-1);
					fprintf(logout,"Age Class\tMales\tMale Deaths\tFathers\tFemales\tFemale Deaths\tMothers\n");
					for(age=0; age<=age_classes; age++) {
													fprintf(logout,"%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\n",age,male_age_dist[age],male_deaths[age],fathers[age],
																					female_age_dist[age],female_deaths[age],mothers[age]);
					}
					fclose(logout);
	/******************************* end of burn in ***************************/
	/* start each repeat from init pop distribution */
					printf("Burn-in age distribution\n"); AgeDist(popN,(long) age_classes,female_age_dist,init_female_age);
					AgeDist(popN,(long) age_classes,male_age_dist,init_male_age);
					printf("Burn-in mut_mortation distribution\n");
					MutDist(popN,(long) age_classes,init_female_age,female_mut_mort_dist,init_female_mut_mort);
					MutDist(popN,(long) age_classes,init_male_age,male_mut_mort_dist,init_male_mut_mort);
					MutDist(popN,(long) age_classes,init_female_age,female_mut_fert_dist,init_female_mut_fert);
					MutDist(popN,(long) age_classes,init_male_age,male_mut_fert_dist,init_male_mut_fert);

					printf("Burn-in loci fixed\n");
					LociFixed(nloci,popN,init_male_mut_mort,male_fixed_mort);
					LociFixed(nloci,popN,init_female_mut_mort,female_fixed_mort);
					LociFixed(nloci,popN,init_male_mut_fert,male_fixed_fert);
					LociFixed(nloci,popN,init_female_mut_fert,female_fixed_fert);

					printf("Burn-in male fertility distribution\n");
					FertDist('m',nloci,popN,age_classes,init_male_age,male_fert_Rdist[0],init_male_mut_fert,male_fertility,age_of_expression,mut_fert_effect);
					printf("Burn-in female fertility distribution\n");
					FertDist('f',nloci,popN,age_classes,init_female_age,female_fert_Rdist[0],init_female_mut_fert,female_fertility,age_of_expression,mut_fert_effect);

					printf("Burn-in male mortality distribution\n");
					MortDist('m',nloci,popN,age_classes,init_male_age,male_mort_Rdist[0],init_male_mut_mort,male_survival,age_of_expression,mut_mort_effect);
					printf("Burn-in female mortality distribution\n");
					MortDist('f',nloci,popN,age_classes,init_female_age,female_mort_Rdist[0],init_female_mut_mort,female_survival,age_of_expression,mut_mort_effect);

					agedistout=fopen(agedistpath,"a");
					if(agedistout == NULL) {printf("Can't open %s\n",agedistpath); exit(1); }
					fprintf(agedistout,"%ld\t",0);
					for(age=0; age<=(age_classes); age++) {fprintf(agedistout,"%ld\t",male_age_dist[age]); }
					for(age=0; age<=(age_classes); age++) {fprintf(agedistout,"%ld\t",female_age_dist[age]); }
					fprintf(agedistout,"\n");
					fclose(agedistout);

					mut_mortdistout=fopen(mut_mortdistpath,"a");
					if(mut_mortdistout == NULL) {printf("Can't open %s\n",mut_mortdistpath); exit(1); }
					fprintf(mut_mortdistout,"%ld\t%ld\t",0,male_fixed_mort[0]);
					for(age=0; age<=(age_classes); age++) {fprintf(mut_mortdistout,"%ld\t",male_mut_mort_dist[age]); }
					fprintf(mut_mortdistout,"%ld\t",female_fixed_mort[0]);
					for(age=0; age<=(age_classes); age++) {fprintf(mut_mortdistout,"%ld\t",female_mut_mort_dist[age]); }
					fprintf(mut_mortdistout,"\n");
					fclose(mut_mortdistout);

					mut_fertdistout=fopen(mut_fertdistpath,"a");
					if(mut_fertdistout == NULL) {printf("Can't open %s\n",mut_fertdistpath); exit(1); }
					fprintf(mut_fertdistout,"%ld\t%ld\t",0,male_fixed_fert[0]);
					for(age=0; age<=(age_classes); age++) {fprintf(mut_fertdistout,"%ld\t",male_mut_fert_dist[age]); }
					fprintf(mut_fertdistout,"%ld\t",female_fixed_fert[0]);
					for(age=0; age<=(age_classes); age++) {fprintf(mut_fertdistout,"%ld\t",female_mut_fert_dist[age]); }
					fprintf(mut_fertdistout,"\n");
					fclose(mut_fertdistout);

					deathdistout=fopen(deathdistpath,"a");
					if(deathdistout == NULL) {printf("Can't open %s\n",deathdistpath); exit(1); }
					fprintf(deathdistout,"%ld\t",0);
					for(age=0; age<=(age_classes); age++) {fprintf(deathdistout,"%ld\t",male_deaths[age]); }
					for(age=0; age<=(age_classes); age++) {fprintf(deathdistout,"%ld\t",female_deaths[age]); }
					fprintf(deathdistout,"\n");
					fclose(deathdistout);

					birthdistout=fopen(birthdistpath,"a");
					if(birthdistout == NULL) {printf("Can't open %s\n",birthdistpath); exit(1); }
					fprintf(birthdistout,"%ld\t",0);
					for(age=0; age<=(age_classes); age++) {fprintf(birthdistout,"%ld\t",fathers[age]); }
					for(age=0; age<=(age_classes); age++) {fprintf(birthdistout,"%ld\t",mothers[age]); }
					fprintf(birthdistout,"\n");
					fclose(birthdistout);

					fertdistout=fopen(fertdistpath,"a");
					if(fertdistout == NULL) {printf("Can't open %s\n",fertdistpath); exit(1); }
					fprintf(fertdistout,"%ld\t",0);
					fprintf(fertdistout,"na\t");
					for(age=1; age<=(age_classes); age++) {
													if(male_fert_Rdist[0][age] >= 0.0) fprintf(fertdistout,"%lf\t",male_fert_Rdist[0][age]);
													else fprintf(fertdistout,"na\t");
					}
					fprintf(fertdistout,"na\t");
					for(age=1; age<=(age_classes); age++) {
													if(female_fert_Rdist[0][age] >= 0.0) fprintf(fertdistout,"%lf\t",female_fert_Rdist[0][age]);
													else fprintf(fertdistout,"na\t");
					}
					fprintf(fertdistout,"\n");
					fclose(fertdistout);

					mortdistout=fopen(mortdistpath,"a");
					if(mortdistout == NULL) {printf("Can't open %s\n",mortdistpath); exit(1); }
					fprintf(mortdistout,"%ld\t",0);
					fprintf(mortdistout,"na\t");
					for(age=1; age<=(age_classes); age++) {
													if(male_mort_Rdist[0][age] >= 0.0) fprintf(mortdistout,"%lf\t",male_mort_Rdist[0][age]);
													else fprintf(mortdistout,"na\t");
					}
					fprintf(mortdistout,"na\t");
					for(age=1; age<=(age_classes); age++) {
													if(female_mort_Rdist[0][age] >= 0.0) fprintf(mortdistout,"%lf\t",female_mort_Rdist[0][age]);
													else fprintf(mortdistout,"na\t");
					}
					fprintf(mortdistout,"\n");
					fclose(mortdistout);

	/*	printf("parents age dist:\n");
	for(age=1;age<=(age_classes-1);age++) {
	for(a=1;a<=(age_classes-1);a++) printf("%ld\t",parents_age_dist[age][a]);
	printf("\n");
	}
	{printf("pause.. ");getchar();printf("\n");}*/
					parentsdistout=fopen(parentspath,"a");
					if(parentsdistout == NULL) {printf("Can't open %s\n",parentspath); exit(1); }
					fprintf(parentsdistout,"Age-Specific distribution of parents' (mother x father) ages\nInitial t=0");
					fprintf(parentsdistout,"Mother's age\t");
					for(age=1; age<=(age_classes-1); age++) fprintf(parentsdistout,"%ld to %ld-\t",(age-1)*age_interval,(age)*age_interval);
					fprintf(parentsdistout,"\n");
					for(age=1; age<=(age_classes-1); age++) {
													fprintf(parentsdistout,"Father's age: %ld to %ld-\t",(age-1)*age_interval,(age)*age_interval);
													for(a=1; a<=(age_classes-1); a++) fprintf(parentsdistout,"%ld\t",parents_age_dist[age][a]);
													fprintf(parentsdistout,"\n");
					}
					fclose(parentsdistout);
					for(age=0; age<=age_classes; age++) {
													for(a=0; a<=age_classes; a++) parents_age_dist[age][a]=0;
					}
	/* parents age distribution taken at end of repeat and averaged over number of repeats*/

	/****** SIMULATION **/
					logout=fopen(logpath,"a");
					if(logout == NULL) {printf("Can't open %s\n",logpath); exit(1); }
					fprintf(logout,"\n********\nSimulation\n");
					fprintf(logout,"Start each repeat with initial age distribution from pre-run\n");
					fprintf(logout,"Simulation time_count\t%ld\n",time_count);
					fprintf(logout,"Simulation total time_count\t%ld\n",total_time_counts);
					fprintf(logout,"Simulation Time (age class interval units)\t%ld\n",time_count*total_time_counts);
					fprintf(logout,"Simulation Repeats\t%ld\n",repeats);
					fclose(logout);
					printf("\nRun %ld repeats with mu_mort = %lf, mut_mort_effect = %lf\n",repeats,mu_mort,mut_mort_effect);
					for(nl=0; nl<=nloci; nl++) {
													male_mut_mort_Ldist[0][nl]=0;
													female_mut_mort_Ldist[0][nl]=0;
													for(i=1; i<=popN; i++) {
																					male_mut_mort_Ldist[0][nl]=male_mut_mort_Ldist[0][nl]+init_male_mut_mort[i][nl];
																					female_mut_mort_Ldist[0][nl]=female_mut_mort_Ldist[0][nl]+init_female_mut_mort[i][nl];
													}
					}
					mystrcpy(mut_mort_locipath,globalpath);
					mystrcat(mut_mort_locipath,run);
					mystrcat(mut_mort_locipath,"/");
					mystrcat(mut_mort_locipath,run);
					mystrcat(mut_mort_locipath,"_locifreq.txt");
					mut_mort_locifile=fopen(mut_mort_locipath,"w");
					if(mut_mort_locifile == NULL) {printf("Can't open %s\n",mut_mort_locipath); exit(1); }
					fprintf(mut_mort_locifile,"Mortality-Causing Mutations Segregating by Locus at end of repeat\n");
					fprintf(mut_mort_locifile,"popN\t%ld\nnloci\t%ld\nmut_mort_rate\t%lf\nmut_mort_effect\t%lf\nRepeats\t%ld\n",popN,nloci,mu_mort,mut_mort_effect,repeats);
					fprintf(mut_mort_locifile,"Repeat (at end)\tTotal (Male), all loci\t");
					for(nl=1; nl<=nloci; nl++) {fprintf(mut_mort_locifile,"(Male) Locus #%ld\t",nl); }
					fprintf(mut_mort_locifile,"Total (Female), all loci\t");
					for(nl=1; nl<=nloci; nl++) {fprintf(mut_mort_locifile,"(Female) Locus #%ld\t",nl); }
					fprintf(mut_mort_locifile,"\n");
					fprintf(mut_mort_locifile,"0\t%ld\t",male_mut_mort_Ldist[0][0]);
					for(nl=1; nl<=nloci; nl++) {fprintf(mut_mort_locifile,"%ld\t",male_mut_mort_Ldist[0][nl]); }
					fprintf(mut_mort_locifile,"%ld\t",female_mut_mort_Ldist[0][0]);
					for(nl=1; nl<=nloci; nl++) {fprintf(mut_mort_locifile,"%ld\t",female_mut_mort_Ldist[0][nl]); }
					fprintf(mut_mort_locifile,"\n");
					fclose(mut_mort_locifile);

	/* begin loop of repeats */
					for(r=1; r<=repeats; r++) {
													printf("\nRepeat #%ld with mu_mort = %lf, mut_mort_effect = %lf\n",r,mu_mort,mut_mort_effect);

	/* start with initial distribution of females */
													for(i=1; i<=popN; i++) {
																					female_age[i] = init_female_age[i];
																					for(nl=0; nl<=nloci; nl++) {female_mut_mort[i][nl] = init_female_mut_mort[i][nl]; }
																					male_age[i] = init_male_age[i];
																					for(nl=0; nl<=nloci; nl++) {male_mut_mort[i][nl] = init_male_mut_mort[i][nl]; }
													}
													LociFixed(nloci,popN,male_mut_mort,male_fixed_mort);
													male_fixed_mort_dynamic[0][r]=(double) male_fixed_mort[0];
													LociFixed(nloci,popN,female_mut_mort,female_fixed_fert);
													female_fixed_mort_dynamic[0][r]=(double) female_fixed_fert[0];
													MutDist(popN,(long) age_classes,male_age,male_mut_mort_Rdist[r],male_mut_mort);
													male_mut_mort_dynamic[0][r]=(double) male_mut_mort_Rdist[r][0];
													MutDist(popN,(long) age_classes,female_age,female_mut_mort_Rdist[r],female_mut_mort);

													female_mut_mort_dynamic[0][r]=(double) female_mut_mort_Rdist[r][0];


	/********************** population with mut_mortality **************************************/
													cum_new_mut_mort[r]=0;
													for(tc=1; tc<=total_time_counts; tc++) {

																					for(ti=1; ti<=time_count; ti++) {
																													t=((tc-1)*time_count) + ti;

	/*intialize for time unit t */
																													for(age=0; age<=age_classes; age++) {
																																					mothers[age]=0;
																																					fathers[age]=0;
																													}
																													new_mut_mort=0;
	/* initial age distributions for time unit t*/
																													AgeDist(popN,(long) age_classes,female_age_dist,female_age);
																													/*printf("\nr=%ld t=%ld: total females = %ld\n",r,t,female_age_dist[0]);*/
																													/*for(age=0;age<=age_classes;age++) {printf("N[%ld]=%ld\n",age,female_age_dist[age]);}*/
																													AgeDist(popN,(long) age_classes,male_age_dist,male_age);

																													/*for(age=0;age<=age_classes;age++) {printf("N[%ld]=%ld\n",age,male_age_dist[age]);}*/

																													MutDist(popN,(long) age_classes,female_age,female_mut_mort_dist,female_mut_mort);

																													/*for(age=0;age<=age_classes;age++) {printf("M[%ld]=%ld\n",age,female_mut_mort_dist[age]);}*/
																													MutDist(popN,(long) age_classes,male_age,male_mut_mort_dist,male_mut_mort);

																													/*for(age=0;age<=age_classes;age++) {printf("M[%ld]=%ld\n",age,male_mut_mort_dist[age]);}*/

	/* determine the number of female deaths in time unit t */
																													Deaths(3,nloci,&idum,age_classes,popN,female_survival,female_mut_mort,age_of_expression,mut_mort_effect,female_age,female_deaths);

																													/*for(age=0;age<=age_classes;age++) {printf("D[%ld]=%ld\n",age,female_deaths[age]);}*/
	/* model 3 is unused */
	/* model will determine the effect of mut_mortalitys on the basis survival table = [male/female]_survival */
	/* determine the number of male deaths in time unit t */

																													Deaths(3,nloci,&idum,age_classes,popN,male_survival,male_mut_mort,age_of_expression,mut_mort_effect,male_age,male_deaths);

																													/*for(age=0;age<=age_classes;age++) {printf("D[%ld]=%ld\n",age,male_deaths[age]);}*/
	/* deaths are indicated by age[i] = 0 */

																													for(i=1; i<=popN; i++) {
																																					if(female_age[i] == 0) {
	/* model 0 => mut_mortalitys affect only females (mothers) */
																																													Parents(2,&idum,popN,nloci,age_classes,i,0,mate,female_age,male_age,female_fertility,male_fertility,
																																																					female_mut_fert,male_mut_fert,age_of_expression,mut_fert_effect,
																																																					&mother,&father);
																																													mothers[female_age[mother]]++;
																																													mothers[0]++;
																																													fathers[male_age[father]]++;
																																													fathers[0]++;
	/* transmission of mut_mortalitys from parents to offspring */
																																													Transmission(nloci,male_mut_mort[father],female_mut_mort[mother],female_mut_mort[i],&idum);
	/* new mut_mortalitys */
																																													new_mut_mort=new_mut_mort+GetMut(&idum,nloci,female_mut_mort[i],mu_mort);
																																					}
																																					if(male_age[i] == 0) {
																																													Parents(2,&idum,popN,nloci,age_classes,0,i,mate,female_age,male_age,female_fertility,male_fertility,
																																																					female_mut_fert,male_mut_fert,age_of_expression,mut_fert_effect,
																																																					&mother,&father);
																																													mothers[female_age[mother]]++;
																																													mothers[0]++;
																																													fathers[male_age[father]]++;
																																													fathers[0]++;
																																													if(t == time_count*total_time_counts) parents_age_dist[male_age[father]][female_age[mother]]++;
	/* transmission of mut_mortalitys from parents to offspring */
																																													Transmission(nloci,male_mut_mort[father],female_mut_mort[mother],male_mut_mort[i],&idum);
	/* new mut_mortalitys */
																																													new_mut_mort=new_mut_mort+GetMut(&idum,nloci,male_mut_mort[i],mu_mort);
																																					}

																													}/*for(i=1;i<=popN;i++) */

	/* update population, deaths have alread been considered, survivors move to next age group */
																													for(i=1; i<=popN; i++) {
																																					if(female_age[i] >= age_classes) {printf("female #%ld age=%ld error1\n",i,female_age[i]); exit(1); }
																																					female_age[i]=female_age[i]+1;
																																					if(male_age[i] >= age_classes) {printf("male #%ld age=%ld error\n",i,male_age[i]); exit(1); }
																																					male_age[i]=male_age[i]+1;
	/*  individual i alread has mut_mortalitys; if i is born in interval, it moves from age 0 to 1, but the offspring mut_mortalitys already belong to i */

																													}/*for(i=1;i<=popN;i++)*/

																					}/*for(ti=1;ti<=time_count;ti++)*/
																					printf("Begin interval: total female mut_mortality = %ld\n",female_mut_mort_dist[0]);
																					printf("End interval: total female mut_mortality = %ld\n",female_mut_mort_Rdist[r][0]);
																					printf("Begin interval: total male mut_mortality = %ld\n",male_mut_mort_dist[0]);
																					printf("End interval: total male mut_mortality = %ld\n",male_mut_mort_Rdist[r][0]);

																					printf("\nr=%ld t=%ld: male deaths = %ld female deaths = %ld total = %ld\n",r,t,male_deaths[0],female_deaths[0],male_deaths[0]+female_deaths[0]);
																					printf("new mut_mortality = %ld\n",new_mut_mort);
																					cum_new_mut_mort[r]=cum_new_mut_mort[r]+new_mut_mort;
																					printf("cumulative new mut_mortality (repeat %ld) = %ld\n",r,cum_new_mut_mort[r]);
																					printf("male demographics repeat %ld tc=%ld\n",r,tc);
																					AgeDist(popN,(long) age_classes,male_age_Rdist[r],male_age);
																					MutDist(popN,(long) age_classes,male_age,male_mut_mort_Rdist[r],male_mut_mort);
																					LociFixed(nloci,popN,male_mut_mort,male_fixed_mort);
																					FertDist('m',nloci,popN,age_classes,male_age,male_fert_Rdist[r],male_mut_fert,male_fertility,age_of_expression,mut_fert_effect);
																					MortDist('m',nloci,popN,age_classes,male_age,male_mort_Rdist[r],male_mut_mort,male_survival,age_of_expression,mut_mort_effect);
																					male_fixed_mort_dynamic[tc][r]=(double) male_fixed_mort[0];
																					male_mut_mort_dynamic[tc][r]=(double) male_mut_mort_Rdist[r][0];

																					printf("female demographics repeat %ld tc=%ld\n",r,tc);
																					AgeDist(popN,(long) age_classes,female_age_Rdist[r],female_age);
																					MutDist(popN,(long) age_classes,female_age,female_mut_mort_Rdist[r],female_mut_mort);
																					LociFixed(nloci,popN,female_mut_mort,female_fixed_mort);
																					FertDist('f',nloci,popN,age_classes,female_age,female_fert_Rdist[r],female_mut_fert,female_fertility,age_of_expression,mut_fert_effect);
																					MortDist('f',nloci,popN,age_classes,female_age,female_mort_Rdist[r],female_mut_mort,female_survival,age_of_expression,mut_mort_effect);
																					female_fixed_mort_dynamic[tc][r]=(double) female_fixed_mort[0];
																					female_mut_mort_dynamic[tc][r]=(double) female_mut_mort_Rdist[r][0];
																					printf("r=%ld tc=%ld female_fixed_mort_dynamic[tc][r]= %lf\n",r,tc,female_fixed_mort_dynamic[tc][r]);
																					/*{printf("pause.. t=%ld ti=%ld",t,ti);getchar();printf("\n");}*/
													}/*** for(tc=1;tc<=total_time_counts;tc++) ***/

	/********* end simulation to t=tmax *********/
	/* loci fixed */
													printf("\nstore data for r=%ld\n",r);
													agedistout=fopen(agedistpath,"a");
													if(agedistout == NULL) {printf("Can't open %s\n",agedistpath); exit(1); }
													fprintf(agedistout,"%ld\t",r);
													for(age=0; age<=(age_classes); age++) {fprintf(agedistout,"%ld\t",male_age_Rdist[r][age]); }
													for(age=0; age<=(age_classes); age++) {fprintf(agedistout,"%ld\t",female_age_Rdist[r][age]); }
													fprintf(agedistout,"\n");
													fclose(agedistout);

													mut_mortdistout=fopen(mut_mortdistpath,"a");
													if(mut_mortdistout == NULL) {printf("Can't open %s\n",mut_mortdistpath); exit(1); }
													fprintf(mut_mortdistout,"%ld\t%ld\t",r,male_fixed_mort[0]);
													for(age=0; age<=(age_classes); age++) {fprintf(mut_mortdistout,"%ld\t",male_mut_mort_Rdist[r][age]); }
													fprintf(mut_mortdistout,"%ld\t",female_fixed_fert[0]);
													for(age=0; age<=(age_classes); age++) {fprintf(mut_mortdistout,"%ld\t",female_mut_mort_Rdist[r][age]); }
													fprintf(mut_mortdistout,"\n");
													fclose(mut_mortdistout);

													deathdistout=fopen(deathdistpath,"a");
													if(deathdistout == NULL) {printf("Can't open %s\n",deathdistpath); exit(1); }
													fprintf(deathdistout,"%ld\t",r);
													for(age=0; age<=(age_classes); age++) {fprintf(deathdistout,"%ld\t",male_deaths[age]); }
													for(age=0; age<=(age_classes); age++) {fprintf(deathdistout,"%ld\t",female_deaths[age]); }
													fprintf(deathdistout,"\n");
													fclose(deathdistout);

													birthdistout=fopen(birthdistpath,"a");
													if(birthdistout == NULL) {printf("Can't open %s\n",birthdistpath); exit(1); }
													fprintf(birthdistout,"%ld\t",r);
													for(age=0; age<=(age_classes); age++) {fprintf(birthdistout,"%ld\t",fathers[age]); }
													for(age=0; age<=(age_classes); age++) {fprintf(birthdistout,"%ld\t",mothers[age]); }
													fprintf(birthdistout,"\n");
													fclose(birthdistout);

													fertdistout=fopen(fertdistpath,"a");
													if(fertdistout == NULL) {printf("Can't open %s\n",fertdistpath); exit(1); }
													fprintf(fertdistout,"%ld\t",r);
													fprintf(fertdistout,"na\t");
													for(age=1; age<=(age_classes); age++) {
																					if(male_fert_Rdist[r][age] >= 0.0) fprintf(fertdistout,"%lf\t",male_fert_Rdist[r][age]);
																					else fprintf(fertdistout,"na\t");
													}
													fprintf(fertdistout,"na\t");
													for(age=1; age<=(age_classes); age++) {
																					if(female_fert_Rdist[r][age] >= 0.0) fprintf(fertdistout,"%lf\t",female_fert_Rdist[r][age]);
																					else fprintf(fertdistout,"na\t");
													}
													fprintf(fertdistout,"\n");
													fclose(fertdistout);

													mortdistout=fopen(mortdistpath,"a");
													if(mortdistout == NULL) {printf("Can't open %s\n",mortdistpath); exit(1); }
													fprintf(mortdistout,"%ld\t",r);
													fprintf(mortdistout,"na\t");
													for(age=1; age<=(age_classes); age++) {
																					if(male_mort_Rdist[r][age] >= 0.0) fprintf(mortdistout,"%lf\t",male_mort_Rdist[r][age]);
																					else fprintf(mortdistout,"na\t");
													}
													fprintf(mortdistout,"na\t");
													for(age=1; age<=(age_classes); age++) {
																					if(female_mort_Rdist[r][age] >= 0.0) fprintf(mortdistout,"%lf\t",female_mort_Rdist[r][age]);
																					else fprintf(mortdistout,"na\t");
													}
													fprintf(mortdistout,"\n");
													fclose(mortdistout);

													for(nl=0; nl<=nloci; nl++) {
																					male_mut_mort_Ldist[r][nl]=0;
																					female_mut_mort_Ldist[r][nl]=0;
																					for(i=1; i<=popN; i++) {
																													male_mut_mort_Ldist[r][nl]=male_mut_mort_Ldist[r][nl]+male_mut_mort[i][nl];
																													female_mut_mort_Ldist[r][nl]=female_mut_mort_Ldist[r][nl]+female_mut_mort[i][nl];
																					}
													}
													mut_mort_locifile=fopen(mut_mort_locipath,"a");
													if(mut_mort_locifile == NULL) {printf("Can't open %s\n",mut_mort_locipath); exit(1); }
													fprintf(mut_mort_locifile,"%ld\t%ld\t",r,male_mut_mort_Ldist[r][0]);
													for(nl=1; nl<=nloci; nl++) {fprintf(mut_mort_locifile,"%ld\t",male_mut_mort_Ldist[r][nl]); }
													fprintf(mut_mort_locifile,"%ld\t",female_mut_mort_Ldist[r][0]);
													for(nl=1; nl<=nloci; nl++) {fprintf(mut_mort_locifile,"%ld\t",female_mut_mort_Ldist[r][nl]); }
													fprintf(mut_mort_locifile,"\n");
													fclose(mut_mort_locifile);

													printf("Finished r=%ld\n",r);
													/*{printf("\nr=%ld pause.. t=%ld ti=%ld",r,t,ti);getchar();printf("\n");}*/
					}/*for(r=1;r<=repeats;r++) */
	/* enter time course for repeat averages */
					/* file for time course (change as you wish) */
					mystrcpy(outpath_time,globalpath);
					mystrcat(outpath_time,run);
					mystrcat(outpath_time,"/");
					mystrcat(outpath_time,run);
					mystrcat(outpath_time,"_time.txt");
					outfile_time=fopen(outpath_time,"w");
					if(outfile_time == NULL) {printf("Can't open %s\n",outpath_time); exit(1); }
					fprintf(outfile_time,"Dynamics\tMeans over repeats\n");
					fprintf(outfile_time,"popN\t%ld\nnloci\t%ld\nmut_mort_rate\t%lf\nmut_mort_effect\t%lf\nRepeats\t%ld\n",popN,nloci,mu_mort,mut_mort_effect,repeats);
					fprintf(outfile_time,"Time\tMean Male mortality loci fixed\tVariance Male mortality loci fixed\t");
					fprintf(outfile_time,"Mean Female mortality loci homozygous in population\tVariance Female mortality loci homozygous in population\t");
					fprintf(outfile_time,"Mean Male mortality loci homozygous in population\tVariance Male mortality loci homozygous in population\t");
					fprintf(outfile_time,"Mean Female mortality mutations segregating in population\tVariance Female mortality mutations segregating in population\n");
					for(tc=0; tc<=total_time_counts; tc++) {
													t=tc*time_count;
													fprintf(outfile_time,"%ld\t",t);
													MeanVar(repeats,male_fixed_mort_dynamic[tc],&mean,&var);
													fprintf(outfile_time,"%lf\t%lf\t",mean,var);
													MeanVar(repeats,female_fixed_mort_dynamic[tc],&mean,&var);
													fprintf(outfile_time,"%lf\t%lf\t",mean,var);
													MeanVar(repeats,male_mut_mort_dynamic[tc],&mean,&var);
													fprintf(outfile_time,"%lf\t%lf\t",mean,var);
													MeanVar(repeats,female_mut_mort_dynamic[tc],&mean,&var);
													fprintf(outfile_time,"%lf\t%lf\n",mean,var);
					}

					fclose(outfile_time);

					parentsdistout=fopen(parentspath,"a");
					if(parentsdistout == NULL) {printf("Can't open %s\n",parentspath); exit(1); }
					fprintf(parentsdistout,"\nAge-Specific distribution of parents' (mother x father) ages\nFinal t=%ld\nAveraged over %ld repeats\n",time_count*total_time_counts,repeats);
					fprintf(parentsdistout,"Mother's age\t");
					for(age=1; age<=(age_classes-1); age++) fprintf(parentsdistout,"%ld to %ld-\t",(age-1)*age_interval,(age)*age_interval);
					fprintf(parentsdistout,"\n");
					for(age=1; age<=(age_classes-1); age++) {
													fprintf(parentsdistout,"Father's age: %ld to %ld-\t",(age-1)*age_interval,(age)*age_interval);
													for(a=1; a<=(age_classes-1); a++) fprintf(parentsdistout,"%6.1lf\t",((double) parents_age_dist[age][a])/((double) repeats));
													fprintf(parentsdistout,"\n");
					}
					fclose(parentsdistout);

					printf("\nFinished!!\n");
					return 0;
	}/*main*/
	/**********    strings  *********************************/
	int mystrcmp( char s1[], char s2[] )
	{
					int i;

					i= 0;
					while( s2[i] != '\0' ) {
													if ( s1[i] != s2[i] ) {
																					return(0);
													}
													i++;
					}
					if ( s1[i] != '\0' ) return(0);
					return(1);
					/* return 1 if true and strings match */
	}/*end mystrcmp*/
	void str_setup(char str0[], char str1[], char str2[], char str3[], char str4[])
	{
					mystrcpy(str0, str1);
					mystrcat(str0, str2);
					mystrcat(str0, str3);
					mystrcat(str0, str4);
					return;
	}
	/*******************/
	void mystrcat( char str0[], char str[] )
	{
					long i, j;

					i = 0;
					if ( str0[0] != '\0' ) {
													while ( str0[i] != '\0' ) {
																					i++;
													}
					}
					j = 0;
					while ( str[j] != '\0' ) {
													str0[i] = str[j];
													i++;
													j++;
					}
					str0[i] = '\0';
					return;
	}
	/******************/
	double mystrtolf( char S[] )
	{
					long i, j, k, l;
					double power, N, F;
					int left[20], right[20], side, sign=0;

					i=0;
					j=0;
					l=0;
					side = 1;
					if ( S[0] == '+' ) {
													sign = 1;
													i++;
					}
					if ( S[0] == '-' ) {
													sign = -1;
													i++;
					}
					if(sign == 0 ) sign = 1;

					/* printf("\nin=> %s\n", S );
					   getchar();*/
					while( S[i] != '\0' )
					{
													if(S[i] == '.' ) {
																					side = 0;
																					i++;
													}
													if( S[i] <= 47 ) return(-1);
													if( S[i] >= 58 ) return(-1);

													if ( side == 1 ) {
																					left[l] = S[i]-48;
																					/*printf("L[%ld] => %d\n", l, left[l] );*/
																					i++;
																					l++;
													}
													else {
																					right[j] = S[i]-48;
																					/*printf("R[%ld] => %d\n", j, right[j] );*/
																					i++;
																					j++;
													}

					}
					/* printf("fini i = %ld", i );
					   getchar();
					   for(k=l-1; k>=0; k-- ) printf("%d", left[k] );
					   printf(".");
					   for(k=0; k<=j-1; k++ ) printf("%d", right[k]);
					   printf("\n");
					   getchar(); */

					N=0;
					power = 1;
					for(k=l-1; k>=0; k-- ) {
													N = N+ power*( (double) left[k] );
													power = 10*power;
													/* printf("N[%ld]=> %lf\n", k, N );
													   getchar(); */
					}
					F=0.1;
					for(k=0; k<=j-1; k++ ) {
													N = N + F*( (double) right[k] );
													F = F*0.1;
													/*  printf("N[%ld]=> %lf\n", k, N );
													   getchar(); */
					}
					N=N* (double) sign;
					return(N);
	}
	/*********/
	void myltos(long N, char istr[])
	{
					int i, k, n;
					char str[120];
					double x, fract, rem;

					x = (double) N;
					i = 0;
					if (x == 0) {
													istr[0]='0';
													istr[1]='\0';
													return;
					}
					if (x < 0) {
													str[0]='-';
													i=1;
													x=-x;
					}
					while (x > 0) {
													x = x/10;
													fract = 10.0005*modf(x, &rem);
													n = (int) (fract);
													/*printf("value = %lf\tfractional = %lf\tinteger = %lf\tdigit = %d\n",
													   x, fract, rem, n);*/
													str[i] = n + 48;
													i++;
													x = rem;
					}
					n = 0;
					istr[i] = '\0';
					i--;
					for (k=i; k>=0; k--) {
													istr[n] = str[k];
													n++;
					}
					return;
	}
	/**************/
	void mydtos(double N, char istr[], int decimals)
	{
					long i, n;
					char nstr[120];
					double x, fraction, number,factor;

					if (N == 0.0) {
													istr[0]='0';
													istr[1]='\0';
													return;
					}

					if (N < 0.0) {istr[0]='-'; istr[1]='\0'; }
					if (N > 0.0) {istr[0]='\0'; N=N+0.0000001; }
					if (N < 0.0) {N = -N; N=N + 0.0000001; }
					factor=0.555555;
					for (i=1; i<=decimals; i++) factor=factor*0.1000000;
					fraction = modf(N, &number);
					fraction=fraction+factor;
					n = (long) number;

					myltos(n,nstr);
					mystrcat(istr,nstr);
					if(decimals <= 0) return;
					mystrcat(istr,".");
					for (i=1; i<=decimals; i++) {
													x=fraction*10.00000000;

													fraction = modf(x, &number);
													n = (long) number;
													myltos(n,nstr);
													mystrcat(istr,nstr);
					}
					return;
	}/**/
	/****/
	long str_len(char str[])
	{
					long i;

					i = 0;
					if (str[i] == '\0') return(0);
					while (str[i] != '\0') i++;
					return(i);
	}/*end str_len*/
	void mystrcpy( char str0[], char str[] )
	{
					int i;
					i = 0;
					while ( str[i] != '\0' ) {
													str0[i] = str[i];
													i++;
					}
					str0[i] = '\0';
					return;
	}/*end mystrcpy*/
	/************************   end strings **********************/
	/***************************  vector lib *********************/
	char **cmatrix(long nrl, long nrh, long ncl, long nch)
	{
					char  **m;
					long i, nrow = nrh-nrl+1, ncol = nch-ncl+1;

					m = (char **) malloc( (size_t) ( (nrow+1)*sizeof(char*) ) );
					if (!m) {printf("memory allocation failure[1] in cmatrix\n"); exit(1); }
					m += 1;
					m -= nrl;
					m[nrl] = (char *) malloc( (size_t) ( (nrow*ncol+1)*sizeof(char) ) );
					if (!m[nrl]) {printf("memory allocation failure[2] in cmatrix"); exit(1); }
					m[nrl] += 1;
					m[nrl] -= ncl;
					for(i=nrl+1; i<=nrh; i++) m[i] = m[i-1] + ncol;
					return(m);
	}
	/****/
	void free_cmatrix(char **m, long nrl, long nrh, long ncl, long nch)
	{
					free((char*) (m[nrl]+ncl-1));
					free((char*) (m+nrl-1));
	}
	/******/
	char *cvector(long nl, long nh)
	{
					char  *v;

					v = (char *) malloc( (size_t) ( (nh-nl+2)*sizeof(char) ) );
					if (!v) {printf("memory allocation failure in cvector"); exit(1); }
					return(v-nl+1);
	}
	/******/
	void free_cvector(char *v, long nl, long nh)
	{
					free((char*) (v+nl-1));
	}
	/******/
	long *lvector(long nl, long nh)
	{
					long  *v;

					v = (long *) malloc( (size_t) ( (nh-nl+2)*sizeof(long) ) );
					if (!v) {printf("memory allocation failure in lvector\n"); exit(1); }
					return(v-nl+1);
	}
	/******/
	void free_lvector(long *v, long nl, long nh)
	{
					free((char*) (v+nl-1));
	}
	/*****/
	double *dvector(long nl, long nh)
	{
					double  *v;
					v = (double *) malloc( (size_t) ( (nh-nl+2)*sizeof(double) ) );

					if (!v) {printf("memory allocation failure in dvector"); exit(1); }
					return(v-nl+1);
	}
	/******/
	void free_dvector(double *v, long nl, long nh)
	{
					/*printf("free_dvector nl = %ld nh = %ld\n",nl,nh);*/
					free((char*) (v+nl-1));
	}
	/********/
	long **lmatrix(long nrl, long nrh, long ncl, long nch)
	{
					long  **m;
					long i, nrow = nrh-nrl+1, ncol = nch-ncl+1;

					m = (long **) malloc( (size_t) ( (nrow+1)*sizeof(long*) ) );
					if (!m) {printf("memory allocation failure[1] in lmatrix"); exit(1); }
					m += 1;
					m -= nrl;
					m[nrl] = (long *) malloc( (size_t) ( (nrow*ncol+1)*sizeof(long) ) );
					if (!m[nrl]) {printf("memory allocation failure[2] in lmatrix\n"); exit(1); }
					m[nrl] += 1;
					m[nrl] -= ncl;
					for(i=nrl+1; i<=nrh; i++) m[i] = m[i-1] + ncol;
					return(m);
	}
	/****/
	void free_lmatrix(long **m, long nrl, long nrh, long ncl, long nch)
	{
					free((char*) (m[nrl]+ncl-1));
					free((char*) (m+nrl-1));
	}
	/********/
	double **dmatrix(long nrl, long nrh, long ncl, long nch)
	{
					double  **m;
					long i, nrow = nrh-nrl+1, ncol = nch-ncl+1;

					m = (double **) malloc( (size_t) ( (nrow+1)*sizeof(double*) ) );
					if (!m) {printf("memory allocation failure[1] in dmatrix"); exit(1); }
					m += 1;
					m -= nrl;
					m[nrl] = (double *) malloc( (size_t) ( (nrow*ncol+1)*sizeof(double) ) );
					if (!m[nrl]) {printf("memory allocation failure[2] in dmatrix"); exit(1); }
					m[nrl] += 1;
					m[nrl] -= ncl;
					for(i=nrl+1; i<=nrh; i++) m[i] = m[i-1] + ncol;
					return(m);
	}
	/****/
	void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
	{
					free((char*) (m[nrl]+ncl-1));
					free((char*) (m+nrl-1));
	}
	/********/
	int *ivector(long nl, long nh)
	{
					int  *v;

					v = (int *) malloc( (size_t) ( (nh-nl+2)*sizeof(int) ) );
					if (!v) {printf("memory allocation failure in ivector\n"); exit(1); }
					return(v-nl+1);
	}
	/******/
	void free_ivector(int *v, long nl, long nh)
	{
					free((char*) (v+nl-1));
	}
	/****/
	unsigned long *ulvector(long nl, long nh)
	{
					unsigned long  *v;

					v = (unsigned long *) malloc( (size_t) ( (nh-nl+2)*sizeof(unsigned long) ) );
					if (!v) {printf("memory allocation failure in ulvector/n"); exit(1); }
					return(v-nl+1);
	}
	/******/
	void free_ulvector(unsigned long *v, long nl, long nh)
	{
					free((char*) (v+nl-1));
	}
	/**************** end  myvectorlib **********************/
	/**************************************/
	long UniformRndInt(long jlo, long jhi, double x)
	{
	/*	j = jlo + (int) ( (double) (jhi-jlo+1) * ran1 (&idum) );  */
					return(jlo + (long) ( ( (double) (jhi-jlo+1) ) * x ) );

	}/* end UniformRndIn */
	/****************/
	double ran1( long *idumPtr )
	{
					int j;
					long k;
					static long iy = 0;
					static long iv[NTAB];
					double temp;
					/* negative to initialize */
					if ( *idumPtr <= 0 || !iy )
					{
													if ( -(*idumPtr) < 1 ) *idumPtr = 1;
													else *idumPtr = -(*idumPtr);
													for ( j = NTAB+7; j >= 0; j-- )
													{
																					k = (*idumPtr)/IQ;
																					*idumPtr = IA*(*idumPtr-k*IQ) - IR*k;
																					if ( *idumPtr < 0 ) *idumPtr += IM3;
																					if ( j < NTAB ) iv[j] = *idumPtr;
													}
													iy = iv[0];
					}
					k = (*idumPtr)/IQ;
					*idumPtr = IA*(*idumPtr-k*IQ) - IR*k;
					if ( *idumPtr < 0 ) *idumPtr += IM3;
					j = iy/XDIV;
					iy = iv[j];
					iv[j] = *idumPtr;
					if ( (temp = AM1*iy) > RNMX ) return RNMX;
					else return temp;
	}    /* end ran1 */
	/*******/
	long  GetMut(long *idumPtr,long nloci,long *mut,double mu_mort)
	{
					long nl,new_mut;
					double x;
					new_mut=0;
					for(nl=1; nl<=nloci; nl++) {
													if(mut[nl] < 0 || mut[nl] > 2) {printf("GetMut: mut[%ld]=%ld\n",nl,mut[nl]); exit(1); }
													if(mut[nl] < 2) {
																					x=ran1(idumPtr);
																					if(x <= mu_mort) {mut[nl] = mut[nl] +1; new_mut++; }
													}
					}
					mut[0]=0;
					for(nl=1; nl<=nloci; nl++) {mut[0]=mut[0]+mut[nl]; }
					return(new_mut);
	}/*GetMut*/
	void  Transmission(long nloci,long *male_mut,long *female_mut,long *offspring_mut,long *idumPtr)
	{
					long nl;
					double x;
	/* []_mut[individual#][i] i =  1,nloci; i=0 is total for individual*/
	/* initailize offspring mutations */
					for(nl=0; nl<=nloci; nl++) {offspring_mut[nl]=0; }
	/* add mutations transmitted from parent loci to offspring_mut */
					for(nl=1; nl<=nloci; nl++) {
													if(male_mut[nl] > 2 || male_mut[nl] < 0 ) {printf("Transmission: male_mut[locus=%ld] = %ld\n",nl,male_mut[nl]); exit(1); }
													if(male_mut[nl] == 2) {offspring_mut[nl] = offspring_mut[nl] +1; }
													if(male_mut[nl] == 1) {
																					x=ran1(idumPtr);
																					if(x <= 0.5) offspring_mut[nl] = offspring_mut[nl] +1;
													}
													if(female_mut[nl] > 2 || female_mut[nl] < 0 ) {printf("Transmission: female_mut[locus=%ld] = %ld\n",nl,female_mut[nl]); exit(1); }
													if(female_mut[nl] == 2) {offspring_mut[nl] = offspring_mut[nl] +1; }
													if(female_mut[nl] == 1) {
																					x=ran1(idumPtr);
																					if(x <= 0.5) offspring_mut[nl] = offspring_mut[nl] +1;
													}
													if(offspring_mut[nl] > 2 || offspring_mut[nl] < 0 ) {printf("Transmission: offspring_mut[locus=%ld] = %ld\n",nl,offspring_mut[nl]); exit(1); }
					}
					for(nl=1; nl<=nloci; nl++) {offspring_mut[0] = offspring_mut[0] + offspring_mut[nl]; }
					/*	if(male_mut[0] > 0 || female_mut[0] > 0) {
					   for(nl=1;nl<=nloci;nl++) {
					   printf("nl=%ld  male parent mut=%ld female parent mut=%ld  offspring mut=%ld\n",nl,male_mut[nl],female_mut[nl],offspring_mut[nl]);}
					   {printf("Transmission: pause..");getchar();printf("\n");}
					   }*/
					return;
	}/*Transmission*/
	void Deaths(int model,long nloci,long *idumPtr,long age_classes,long popN,double *survival,
									long **mut_mort,long *age_of_expression,double mut_mort_effect,long *age,long *deaths)
	{
					long i;
					double x,surv;
	/* initialize death count */
					for(i=0; i<=age_classes; i++) deaths[i]=0;
	/* determine individuals dying */
					for(i=1; i<=popN; i++) {
													if(age[i] <= 0 || age[i] > age_classes) {printf("age[%ld] = %ld Terminate!\n",i,age[i]); exit(1); }
													/*printf("i=%ld age = %ld\n",i,age[i]);*/
													surv=GetSurvival(nloci,age_classes,age_of_expression,mut_mort_effect,age[i],survival,mut_mort[i]);
													x=ran1(idumPtr);
													/*printf("x=%lf survival[#%ld][age=%ld] = %lf\n",x,i,age[i],surv);*/
													if(x > surv) {deaths[0]++; deaths[age[i]]++; age[i]=0; }
													/*{printf("Deaths: pause..");getchar();printf("\n");}*/
					}
					return;
	}/*Deaths*/
	double GetSurvival(long nloci,long age_classes,long *age_of_expression,double mut_mort_effect,
																long age,double *survival,long *mut_mort)
	{
					long l;
					double surv,d;
					surv=survival[age];
					if(mut_mort[0] <= 0) return(surv);
					for(l=1; l<=nloci; l++) {
													if(age >= age_of_expression[l]) {
																					d=((double) (age-age_of_expression[l]))*mut_mort_effect;
																					if(d > 1.0) d=1.0;
																					surv=(surv)*pow(1.0-d,mut_mort[l]);
													}
					}
					return(surv);
	}/*GetSurvival*/
	double GetFertility(long nloci,long age_classes,long *age_of_expression,double mut_fert_effect,
																	long age,double *fertility,long *mut_fert)
	{
					long l;
					double fert,d;
					fert=fertility[age];
					if(mut_fert[0] <= 0) return(fert);
					for(l=1; l<=nloci; l++) {
													if(age >= age_of_expression[l]) {
																					d=((double) (age-age_of_expression[l]))*mut_fert_effect;
																					if(d > 1.0) d=1.0;
																					fert=(fert)*pow(1.0-d,mut_fert[l]);
													}
					}
					return(fert);
	}/*GetFertility*/
	void Parents(int model,long *idumPtr,long popN,long nloci,long age_classes,long female,long male,double **mate,
										long *female_age, long *male_age, double *female_fertility, double *male_fertility,
										long **female_mut_fert,long **male_mut_fert,long *age_of_expression, double mut_fert_effect,
										long *motherPtr, long *fatherPtr)
	{
					int done,done1;
					long father,mother;
					double x,mother_fert,father_fert;
	/* mode:l (0 female only, 1 male only, 2 both) mutations affect amle & female fertility equally*/
	/*mate=dmatrix(0,age_classes,0,age_classes); mate[male age][female age] */
					if(model != 2) {printf("Parents model must be 2\n"); exit(1); }
					done1=0;
					while(done1 == 0) {
	/* pick female parent at random */
													done=0;
													while(done == 0) {
																					/* pick female parent at random */
																					mother=UniformRndInt(1,popN,ran1(idumPtr));
																					/*printf("pick mother = %ld [female child # %ld]\n",mother,female);*/
																					if(mother != female) {
																													/*printf("mother=%ld age=%ld fertility=%lf\n",mother,female_age[mother],female_fertility[female_age[mother]]);*/
																													mother_fert=GetFertility(nloci,age_classes,age_of_expression,mut_fert_effect,female_age[mother],female_fertility,female_mut_fert[mother]);
																													x=ran1(idumPtr);
																													if(x < mother_fert) {done=1; }
																					}
													}
													/*{printf("pause..");getchar();printf("\n");}*/
	/* pick male parent at random */
													done=0;
													while(done == 0) {
																					/* pick male parent at random */
																					father=UniformRndInt(1,popN,ran1(idumPtr));
																					/*printf("pick father = %ld [male child # %ld]\n",father,male);*/
																					if(father != male) {
																													father_fert=GetFertility(nloci,age_classes,age_of_expression,mut_fert_effect,male_age[father],male_fertility,male_mut_fert[father]);
																													x=ran1(idumPtr);
																													if(x < father_fert) {done=1; }
																					}
													}
													/*{printf("pause..");getchar();printf("\n");}*/
													x=ran1(idumPtr);
													if(x <  mate[male_age[father]][female_age[mother]]) {done1=1; }
													/*{printf("pause done1=%d ..",done1);getchar();printf("\n");}*/
					}/*while(done1 == 0) */
	/*	printf("child male=%ld female=%ld\n",male,female);
	printf("father=%ld mother=%ld\n",father,mother);
	printf("father age=%ld; fertility = %lf\n",male_age[father],father_fert);
	printf("mother age=%ld; fertility = %lf\n",female_age[mother],mother_fert);
	printf("mate[father][mother] = %lf\n",mate[male_age[father]][female_age[mother]]);
	*/
					(*motherPtr) = mother;
					(*fatherPtr) = father;
					return;
	}/*Parents*/
	void AgeDist(long popN, long age_classes,long *age_dist, long *age)
	{
					long i,a;
					for(a=0; a<=age_classes; a++) age_dist[a]=0;
					for(i=1; i<=popN; i++) {
													if(age[i] <= 0 || age[i] > age_classes) {printf("AgeDistribution age[%ld]=%ld error\n",i,age[i]); exit(1); }
													age_dist[0]++;
													age_dist[age[i]]++;
					}
					return;
	}
	void MutDist(long popN, long age_classes, long *age,long *mut_dist, long **mut)
	{
					long i,a;
					for(a=0; a<=age_classes; a++) mut_dist[a]=0;
					for(i=1; i<=popN; i++) {
													if(mut[i][0] < 0) {printf("MutDistribution mut[%ld]=%ld error\n",i,mut[i][0]); exit(1); }
													if(age[i] <= 0 || age[i] > age_classes) {printf("MutDistribution age[%ld]=%ld error\n",i,age[i]); exit(1); }
													mut_dist[0]=mut_dist[0]+mut[i][0];
													/*printf("MutDist age[#%ld] = %ld mut=%ld\n",i,age[i],mut[i][0]);*/
													mut_dist[age[i]]=mut_dist[age[i]]+mut[i][0];
					}
					return;
	}
	void FertDist(char sex,long nloci,long popN, long age_classes, long *age,double *fert_dist,
											long **mut,double *fertility,long *age_of_expression, double mut_effect)
	{
					long i,a;
					double fert,*age_dist;
					age_dist=dvector(0,age_classes);
					for(a=0; a<=age_classes; a++) {fert_dist[a]=0.0,age_dist[a]=0.0; }
					for(i=1; i<=popN; i++) {
													if(age[i] <= 0 || age[i] > age_classes) {printf("FertDist age[%ld]=%ld error\n",i,age[i]); exit(1); }
													/*printf("GetFert %ld ag=%ld mut=%ld\n",i,age[i],mut[i][0]);*/
													fert=GetFertility(nloci,(long)age_classes,age_of_expression,mut_effect,age[i],fertility,mut[i]);
													/*if(mut[i][0] == 2) printf("#%ld [%ld] mut=%ld fertility=%lf fert=%lf\n",i,age[i],mut[i][0],fertility[age[i]],fert);
													   {printf("pause.");getchar();printf("\n");}*/
													fert_dist[age[i]]=fert_dist[age[i]]+fert;
													age_dist[age[i]]++;
					}
					for(a=1; a<=age_classes; a++) {
													if(age_dist[a] > 0.0) {
																					/*printf("a=%ld total fert_dist[a]=%lf ",a,fert_dist[a]);*/
																					fert_dist[a]=fert_dist[a]/age_dist[a];
													} else fert_dist[a] = -1.0;
					}
					free_dvector(age_dist,0,age_classes);
					/*printf("return\n");*/
					return;
	}/*FertDist*/
	void MortDist(char sex,long nloci,long popN, long age_classes, long *age,double *mort_dist,
											long **mut,double *survival,long *age_of_expression, double mut_effect)
	{
					long i,a;
					double mort,*age_dist;
					age_dist=dvector(0,age_classes);
					for(a=0; a<=age_classes; a++) {mort_dist[a]=0.0,age_dist[a]=0.0; }
					for(i=1; i<=popN; i++) {
													if(age[i] <= 0 || age[i] > age_classes) {printf("FertDist age[%ld]=%ld error\n",i,age[i]); exit(1); }
													/*printf("GetFert %ld ag=%ld mut=%ld\n",i,age[i],mut[i][0]);*/
													mort=GetSurvival(nloci,(long)age_classes,age_of_expression,mut_effect,age[i],survival,mut[i]);
													/*{printf("pause.");getchar();printf("\n");}*/
													mort_dist[age[i]]=mort_dist[age[i]]+mort;
													age_dist[age[i]]++;
					}
					for(a=1; a<=age_classes; a++) {
													if(age_dist[a] > 0.0) {
																					/*printf("a=%ld total fert_dist[a]=%lf ",a,fert_dist[a]);*/
																					mort_dist[a]=mort_dist[a]/age_dist[a];
													} else mort_dist[a] = -1.0;
					}
					free_dvector(age_dist,0,age_classes);
					/*printf("return\n");*/
					return;
	}/*MortDist*/
	void  LociFixed(long nloci,long popN,long **mut,long *fixed)
	{
					long i,nl;
					for(nl=0; nl<=nloci; nl++) fixed[nl]=0;
					for(i=1; i<=popN; i++) {
													for(nl=1; nl<=nloci; nl++) {
																					if(mut[i][nl] < 0 || mut[i][nl] > 2) {printf("LociFixed mut[%ld][%ld] = %ld\n",i,nl, mut[i][nl]); exit(1); }
																					if(mut[i][nl] == 2) fixed[nl]++;
													}
					}
					for(nl=1; nl<=nloci; nl++) fixed[0]=fixed[0]+fixed[nl];
					return;
	}/*LociFixed*/
	void MeanVar(long n, double *data, double *meanPtr, double *svarPtr)
	{
					long i;
					double sumx,sumxsqr,mean,f,var,count;

					sumx=0.0;
					count=(double) n;
					f=(count)/(count-1.0);
					sumxsqr=0.0;
					for(i=1; i<=n; i++) {
													sumx=sumx + data[i];
													sumxsqr=sumxsqr + data[i]*data[i];
					}
					mean=sumx/count;
					var=(sumxsqr/count) - (mean*mean);
					*meanPtr=mean;
					*svarPtr=f*var;
					return;

	}/* end MeanVar */
