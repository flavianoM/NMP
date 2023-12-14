/* Copyright: Flaviano Morone */ 
/* Created on July 6th, 2023 */

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<time.h>

#define N 1024
#define TOL 0.001
#define T1 500
#define T2 500

char fname[N];
int A[N][N], deg_row[N], deg_col[N], A_row[N][N], A_col[N][N];
double r[N], s[N], z[N], w[N], v[N], nu[N];
double r_new[N], s_new[N];
double v_new[N], nu_new[N];
int P[N], Q[N], R[N], S[N];
int p[N][N], q[N][N];


//import bipartite matrix "A"
void import_network(char *filename){
    int n,m,aij,i,a;
    FILE *fin;
    char line[1024], *token;
    double norm;

    sprintf(fname, "%s", filename);
    fin = fopen(fname, "r");
    if(fin == NULL){
        perror("Error loading network");
        exit(0);
    }

    n=0;
    while(fgets(line, 1024, fin) != NULL){
        n++;
        token = strtok(line,",");
        m=0;
        while(token != NULL){
            m++;
            aij = atoi(token);
            A[n][m] = aij > 0 ? 1 : 0;   //binarize
            token = strtok(NULL, ",");
        }
    }
    fclose(fin);
    A[0][0]=n;
    A[0][1]=m;

    printf("%s imported.\n N (rows) = %d\n M (columns) = %d\n", filename, A[0][0], A[0][1]);
}

//create adjacency lists "Ain", "Aout"
void make_network(void){
    int n,m;
    int i,ni,a,na;
    n=A[0][0];
    m=A[0][1];

    for(i=1;i<=n;i++){
        deg_row[i]=0;
        for(a=1;a<=m;a++){
            deg_row[i] += A[i][a];
        }
    }
    for(a=1;a<=m;a++){
        deg_col[a]=0;
        for(i=1;i<=n;i++){
            deg_col[a] += A[i][a];
        }
    }
    for(i=1;i<=n;i++){
        na=1;
        for(a=1;a<=m;a++){
            if(A[i][a]>0){
                A_row[i][na]=a;
                na++;
            }
        }
    }
    for(a=1;a<=m;a++){
        ni=1;
        for(i=1;i<=n;i++){
            if(A[i][a]>0){
                A_col[a][ni]=i;
                ni++;
            }
        }
    }
}

//initialize ranking vectors "r", "s"
//initialize sinkhorn vectors "v", "nu"
void init_all(int n, int m, double beta) {
    int i, a, ni, na;
    double sum;

    for(i = 1; i <= n; i++) {
        r[i] = i;
        v[i] = drand48();
    }
    for(a = 1; a <= m; a++) {
        s[a] = a;
        nu[a] = drand48();
    }
    //init z[i] and w[a]
    for(i = 1; i <= n; i++) {
        sum = 0;
        for(a = 1; a <= deg_row[i]; a++) {
            na = A_row[i][a];
            sum += s[na];
        }
        z[i] = exp(-beta*sum);
    } 
    for(a = 1; a <= m; a++) {
        sum = 0;
        for(i = 1; i <= deg_col[a]; i++) {
            ni = A_col[a][i];
            sum += r[ni];
        }
        w[a] = exp(-beta*sum);
    }
}

//sub-routine to update sinkhorn vectors
double sinkhorn(int n, int m){
    int i, j, k, t;
    int a, b, c;
    double zi, wa;
    double sum, SUM, delta, diff;

    t = 1;
    do {
        diff = 0.0;
        for(j = 1; j <= n; j++) {
            SUM = 0;
            for(i = 1; i <= n; i++) {
                sum = 0;
                zi = 1.0/pow(z[i], j);
                for(k = 1; k <= n; k++) {
                    zi *= z[i];
                    sum += zi * v[k];
                }
                SUM += 1.0/sum;
            }
            v_new[j] = 1.0/SUM;
            delta = fabs(v_new[j] - v[j]);
            diff = (delta > diff ? delta : diff);
        }
        
        for(b = 1; b <= m; b++) {
            SUM = 0;
            for(a = 1; a <= m; a++) {
                sum = 0;
                wa = 1.0/pow(w[a], b);
                for(c = 1; c <= m; c++) {
                    wa *= w[a];
                    sum += wa * nu[c];
                }
                SUM += 1.0/sum;
            }
            nu_new[b] = 1.0/SUM;
            delta = fabs(nu_new[b] - nu[b]);
            diff = (delta > diff ? delta : diff);
        }
        for(j = 1; j <= n; j++){
            v[j] = v_new[j];
        }
        for(b = 1; b <= m; b++) {
            nu[b] = nu_new[b];
        }
        t++;
    }while(diff > TOL && t < T1);
    
    return diff;
}

//one iteration of the NMP
double onestep(double beta, int n, int m, double eps) {
    int i, j, a, b, na, ni; 
    double zi, wa;
    double sum, num, den, x, delta, diff;

   // update z[i] and w[a]
    for(i = 1; i <= n; i++) {
        sum = 0;
        for(a = 1; a <= deg_row[i]; a++) {
            na = A_row[i][a];
            sum += s[na];
        }
        z[i] = exp(-beta*sum);
    } 
    for(a = 1; a <= m; a++) {
        sum = 0;
        for(i = 1; i <= deg_col[a]; i++) {
            ni = A_col[a][i];
            sum += r[ni];
        }
        w[a] = exp(-beta*sum);
    }
    // sinkhorn step
    sinkhorn(n, m);

    diff = 0;
    for(i = 1; i <= n; i++){
        num = 0;
        den = 0;
        zi = 1;
        for(j = 1; j <= n; j++){
            zi *= z[i];
            x = zi * v[j];
            num += j*x;
            den += x;
        }
        r_new[i] = eps*(num/den) + (1-eps)*r[i];
        delta = fabs(r_new[i] - r[i]);
        diff = (delta > diff ? delta : diff);
        r[i] = r_new[i];
    }
    for(a = 1; a <= m; a++) {
        num = 0; 
        den = 0;
        wa = 1;
        for(b = 1; b <= m; b++) {
            wa *= w[a];
            x = wa * nu[b];
            num += b*x;
            den += x; 
        }
        s_new[a] = eps*(num/den) + (1-eps)*s[a];
        delta = fabs(s_new[a] - s[a]);
        diff = (delta > diff ? delta : diff);
        s[a] = s_new[a];
    }
    return diff;
}
/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
//////////////////// Mergesort (I need it to sort the rankings)//////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
typedef struct{
	double fit;
	int id;
}FIT;

//mergesort
void merge(FIT *a, FIT *b, FIT *c, int m, int n) {
	int i, j, k;
	i = 0;
	j = 0;
	k = 0;
    while(i < m && j < n) {
		if(a[i].fit < b[j].fit) { 
			c[k].fit = a[i].fit;
			c[k].id = a[i].id;
			k++;
			i++;
		}
		else {
			c[k].fit = b[j].fit;
			c[k].id = b[j].id;
			k++;
			j++;
		}
	}
	while(i<m) {
		c[k].fit = a[i].fit;
		c[k].id = a[i].id;
		k++;
		i++;
	}
    while(j<n) {
		c[k].fit = b[j].fit;
		c[k].id = b[j].id;
		k++;
		j++;
	}
}
void sort(FIT *F, int n) {
	int j, k, a = 1, inc = 0, cnt=0, lenght = 0, temp = n;
    FIT *wf, *yf;
    while(temp != 0) {
		while(a < temp - a)
			a *=2;
		lenght +=a;
		wf = calloc(a, sizeof(FIT));		
		for(k = 1; k < a; k *= 2){
			for(j = 0; j < a - k; j += 2 * k) {
				merge(F + j + inc, F + j + k + inc, wf + j, k, k);
			}
			for(j = 0; j < a; ++j) {
				F[j + inc].fit = wf[j].fit;
				F[j + inc].id = wf[j].id;
			}
		}
		++cnt;
		free(wf);
		inc +=a;
		temp -=a;	   
		if(cnt >= 2 ) {
			yf = calloc(lenght, sizeof(FIT));
			merge(F, F + lenght-a, yf, lenght-a, a);
			for(j = 0; j < lenght; ++j) {
				F[j].fit = yf[j].fit;
				F[j].id = yf[j].id;
			}			
			free(yf);
		}
		a = 1;
    }
}
//////////////////// End of Mergesort ///////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////


//compute energy
double get_energy(int n, int m){
    int i, a;
    double sum, norm;

    sum = 0;
    norm = 0;    
    for(i = 1; i <= n; i++) {
        for(a = 1; a <= m; a++) {
            sum += R[i]*S[a]*A[i][a];
            norm += A[i][a];
        }
    }
    sum = sum/(n*m*norm);
    return sum;
}

inline int maxint(int a1, int a2){
    return (a1>a2?a1:a2);
}
inline double maxfloat(double a1, double a2){
    return (a1>a2?a1:a2);
}
int get_ki_max(int n){
    int i, ki;
    ki = 0;
    for(i = 1; i <= n; i++) {
        ki = (deg_row[i] > ki ? deg_row[i] : ki);
    }
    return ki;
}
int get_ka_max(int m){
    int a, ka;
    ka = 0;
    for(a = 1; a <= m; a++) {
        ka = (deg_col[a] > ka ? deg_col[a] : ka);
    }
    return ka;
}

//get permutations "P" and "Q" to rearrange adj matrix
void get_P(int n) {
    int i, k;
    
    k = 1;
    while(k <= n) {
        for(i = 1; i <= n; i++) {
            if(R[i] == k){
                P[k] = i;
            }
        }
        k++;
    }
}
void get_Q(int m) {
    int a, c;
    
    c = 1;
    while(c <= m) {
        for(a = 1; a <= m; a++) {
            if(S[a] == c){
                Q[c] = a;
            }
        }
        c++;
    }
}


/////////////////////////////////////////// Main Program ///////////////////////////////////////////////////////
int main(int argc, char *argv[]){
    int n,m;
    int i,j,k,a,b,t,ki,ka;
    double beta, diff, eps, bmax, beta1, e_gs;
    char *string;
    FIT *fit;
    FILE *fout;

    string = argv[1];
    beta1 = 20.0;    //choose bigger value for better result (convergence may slow down)
    srand48(time(NULL));

    import_network(string);

     make_network();
    n = A[0][0];
    m = A[0][1];
    fit = (FIT *)calloc(N+1, sizeof(FIT));

    ki = get_ki_max(n);
    ka = get_ka_max(m);
    bmax = maxfloat(n*ki, m*ka);

    beta = beta1/bmax;
    eps = 1.0;
    init_all(n, m, beta);

    t = 1;
    do {
        diff = onestep(beta, n, m, eps);
        t++;
    } while(diff > TOL && t < T2);
    
    for(i = 0; i < n; i++) {
        fit[i].fit = r[i+1];
        fit[i].id = i+1;
    }
    sort(fit, n);

    for(i = 0; i < n; i++) {
        R[fit[i].id] = i+1;
    }
    for(a = 0; a < m; a++) {
        fit[a].fit = s[a+1];
        fit[a].id = a+1;
    }
    sort(fit, m);

    for(a = 0; a < m; a++) {
        S[fit[a].id] = a+1;
    }

    //print ground state energy (minimal cost)
    e_gs = get_energy(n, m);
    printf("E_gs = %f\n", e_gs);
    fflush(stdout);

    get_P(n);
    get_Q(m);


    //print nested matrix 
    sprintf(fname, "Anest_%s.txt", string );
    fout = fopen(fname, "w");
    for(i = 1; i <= n; i++){
        for(a = 1; a <= m; a++){
            if(A[P[i]][Q[a]] > 0) {
                fprintf(fout, "%d %d %d\n", a, n-i+1, A[P[i]][Q[a]]);
            }
        }
    }
    fclose(fout);


    //print a summary of the results
    sprintf(fname, "result_%s.txt", string);
    fout = fopen(fname, "w");
    fprintf(fout, "beta = %f\n", beta);
    fprintf(fout, "ROWS = %d\n", n);
    fprintf(fout, "COLS = %d\n", m);
    fprintf(fout, "E_gs = %f\n", e_gs);
    fprintf(fout, "\n");
    fprintf(fout, "ROW RANK\n");
    fprintf(fout, "i, R[i]\n");
    for(i = 1; i <= n; i++){
        fprintf(fout, "%d %d\n", i , R[i]);
    }
    fprintf(fout, "\n");

    fprintf(fout, "COL RANK\n");
    fprintf(fout, "a, S[a]\n");
    for(a = 1; a <= m; a++){
        fprintf(fout, "%d %d\n", a , R[a]);
    }
    fclose(fout);









    return 0;
}
