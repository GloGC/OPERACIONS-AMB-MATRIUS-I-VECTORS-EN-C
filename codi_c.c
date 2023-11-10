#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#define N 512

float Mat[N][N];
float MatDD[N][N];
float V1[N];
float V2[N];
float V3[N];
float V4[N];
void InitData() {
        int i, j;
        srand(8824553);
        for (i = 0; i < N; i++) {
                for (j = 0; j < N; j++) {
                        Mat[i][j] = (((i * j) % 3) ? -1 : 1) * (100.0 * (rand() / (1.0 * RAND_MAX)));
                        if ((abs(i - j) <= 3) && (i != j))
                                MatDD[i][j] = (((i * j) % 3) ? -1 : 1) * (rand() / (1.0 * RAND_MAX));
                        else if (i == j)
                                MatDD[i][j] = (((i * j) % 3) ? -1 : 1) * (10000.0 * (rand() / (1.0 * RAND_MAX)));
                        else
                                MatDD[i][j] = 0.0;
                }
        V1[i] = (i < N / 2) ? (((i * j) % 3) ? -1 : 1) * (100.0 * (rand() / (1.0 * RAND_MAX))) : 0.0;
        V2[i] = (i >= N / 2) ? (((i * j) % 3) ? -1 : 1) * (100.0 * (rand() / (1.0 * RAND_MAX))) : 0.0;
        V3[i] = (((i * j) % 5) ? -1 : 1) * (100.0 * (rand() / (1.0 * RAND_MAX)));
        }
}

void PrintVect(float vect[N], int from, int numel){
        for (int i = from; i < from + numel; i++)
                printf("%f ", vect[i]);
        printf("\n");
}

void PrintRow(float mat[N][N], int row, int from, int numel){
        for (int i = from; i < from + numel; i++)
                printf("%f ", Mat[row][i]);
        printf("\n");
}

void MultEscalar(float vect[N], float vectres[N], float alfa){
        for (int i = 0; i < N; i++){
                vectres[i] = alfa * vect[i];
                printf("%f", vectres[i]);
        }
        printf(" \n");
}

float Scalar(float V1[N], float V2[N]){
        float prodesc = 0.0;
        for (int i = 0; i < N; i++){
                prodesc += V1[i] * V2[i];
        }

        return prodesc;
}

float Magnitude(float vect[N]){
        float magnitud = 0.0;
        float quadrat = 0.0;
        float sumatori = 0.0;
        for (int i = 0; i < N; i++){
                quadrat = vect[i] * vect[i];
                sumatori += quadrat;
        }
        magnitud = sqrt(sumatori);
        printf("%f", magnitud);
        printf("\n");
}

int Ortogonal(float vect1[N], float vect2[N]){
        float prodesc = 0.0;
        for (int i = 0; i < N; i++){
                prodesc = vect1[i] * vect2[i];
        }
        return prodesc;
        if (prodesc == 0 ){
                return 1;
        }else{
                return 0;
        }
}

void Projection(float vect1[N], float vect2[N], float vectres[N]){
        float pe = Scalar(V1, V2);
        float mg = Magnitude(V2);
        float proj = 0.0;
        for (int i = 0; i < N; i++){
                proj = (pe * V2[i]) / mg;
                printf("%f", proj);
        }
}

float Infininorm(float M[N][N]){
        float sumMax = 0.0;
        float sumfila = 0.0;
        for (int i = 0; i < N; i++){
                        for (int j = 0; j < N; j++){
                                sumfila += fabs(M[i][j]);
                        }
                        if (sumfila > sumMax){
                                sumMax = sumfila;
                        }
        }
        return sumMax;
}

float Onenorm(float M[N][N]){
        float sumMax = 0.0;
        float sumcolumn = 0.0;
        for (int j = 0; j < N; j++){
                        for (int i = 0; i < N; i++){
                                sumcolumn += fabs(M[i][j]);
                        }
                        if (sumcolumn > sumMax){
                                sumMax = sumcolumn;
                        }
        }
        return sumMax;
}

float NormFrobenius(float M[N][N]){
        float Frob = 0.0;
        float sumquadrat = 0.0;
        for (int i = 0; i < N; i++){
                for (int j = 0; j < N; j++){
                        sumquadrat += (M[i][j] * M[i][j]);
                }
                Frob = sqrt(sumquadrat);
        }
        return Frob;

}

int DiagonalDom(float M[N][N]){
        float vad = 0.0;
        float diagonal = 0.0;
        float conterrors = 0.0;
        float sum = 0.0;
        float sumrest = 0.0;
        for (int i = 0; i < N; i++){
                diagonal = M[i][i];
                vad = fabs(diagonal);
                for (int j = 0; j < N; j++){
                        sum += M[i][j];
                        sumrest = (sum - diagonal);
                        if (vad < sumrest){
                                conterrors += 1;
                        }
                }
                if (conterrors == 0){
                        return 1; // La matriu és dominant
                }else{
                        return 0; // La matriu no és dominant
                }
        }
}

int Jacobi(float M[N][N] , float vect[N], float vectres[N], unsigned iter){
        int resultat = DiagonalDom(M);
        if (resultat == 1){
                float x[N];
                for (unsigned k = 0; k < iter; k++){
                        for (int i = 0; i < N; i++){
                                for (int j = 0; j < N; j++){
                                        vectres[i] = ((vect[i] - (Mat[i][j] * x[j] - Mat[i][i])) / Mat[i][i]);
                                }
                        x[i] = vectres[i];
                }
                return 1; // Es pot aplicar el mètode Jacobi
        }
        }else{
                return 0; // No es pot aplicar el mètode Jacobi
        }
}

int main(){
        InitData();

// A)
        printf("V1 del 0 al 9 és");
        PrintVect(V1, 0, 10);
        printf("V1 del 256 al 265 és");
        PrintVect(V1, 256, 10);
        printf("V2 del 0 al 9 és");
        PrintVect(V2, 0, 10);
        printf("V2 del 256 al 265 és");
        PrintVect(V2, 256, 10);
        printf("V3 del 0 al 9 és");
        PrintVect(V3, 0, 10);
        printf("V3 del 256 al 265 és");
        PrintVect(V3, 256, 10);

// B)
        printf("Els elements del 0 al 9 de la fila 0 són");
        PrintRow(Mat, 0, 0, 10);
        printf("Els elements del 0 al 9 de la fila 100 són");
        PrintRow(Mat, 100, 0, 10);

// c)
        printf("Els elements del 0 al 9 de la fila 0");
        PrintRow(MatDD, 0, 0, 10);
        printf("Els elements del 90 al 99 de la fila 100 són");
        PrintRow(MatDD, 100, 90, 10);

// H)
        printf("El resultat del producte de V3 amb l'escalar 2 del 0 al 9 és");
        MultEscalar(V3, V3, 2);
        PrintVect(V3, 0, 10);
        printf("El resultat del producte de V3 amb l'escalar 2 del 256 al 265 és");
        MultEscalar(V3, V3, 2);
        PrintVect(V3, 256, 10);

// E)
        printf("El producte escalar de V1 i V2 és %f.\n", Scalar(V1, V2));
        printf("El producte escalar de V1 i V2 és %f.\n", Scalar(V1, V3));
        printf("El producte escalar de V2 i V3 és %f.\n", Scalar(V2, V3));

// F)
        printf("La magintud de V1");
        Magnitude(V1);
        printf("La magintud de V2 és");
        Magnitude(V2);
        printf("La magintud de V3 és");
        Magnitude(V3);

// G)
        int x = Ortogonal(V1, V2);
        if (x == 1){
                printf("els vectors V1 i V2 són ortogonals");
        }else{
                printf("els vectors V1 i V2 no són ortogonals");
        }
        int y = Ortogonal(V1, V3);
        if (y == 1){
                printf("els vectors V1 i V3 són ortogonals");
        }else{
                printf("els vectors V1 i V3 no són ortogonals");
        }
        int z = Ortogonal(V2, V3);
        if (z == 1){
                printf("els vectors V2 i V3 són ortogonals");
        }else{
                printf("els vectors V2 i V3 no són ortogonals");
        }

// I)
        printf("Els 10 primers elements de la projecció de V2 sobre V3 són");
        Projection(V2, V3, V4);
        printf("Els 10 primers elements de la projecció de V1 sobre V2 són");
        Projection(V1, V2, V4);

// D)
        printf("La infininorma de Mat és %f.\n", Infininorm(Mat));
        printf("La norma-ú de Mat és %f.\n", Onenorm(Mat));
        printf("La Norma de Frobenius de Mat és %f.\n", NormFrobenius(Mat));
        int Diagx = DiagonalDom(Mat);
        if (Diagx == 1){
                printf("La matriu Matt és diagonal dominant");
        }else{
                printf("La matriu Matt no és diagonal dominant");
        }

        printf("La infininorma de Mat és %f.\n", Infininorm(MatDD));
        printf("La norma-ú de Mat és %f.\n", Onenorm(MatDD));
        printf("La Norma de Frobenius de Mat és %f.\n", NormFrobenius(MatDD));
        int Diagy = DiagonalDom(MatDD);
        if (Diagy == 1){
                printf("La matriu MattDD és diagonal dominant");
        }else{
                printf("La matriu MattDD no és diagonal dominant");
        }

// J)
        printf("Els deu primers elements de la solució aproximada feta amb 1 iteració de MatDD són");
        Jacobi(MatDD, V3, V4, 1 );
        PrintVect(V4, 0, 10);
        printf("Els deu primers elements de la solució aproximada feta amb 1000 iteracions de MatDD són");
        Jacobi(MatDD, V3, V4, 1000);
        PrintVect(V4, 0, 10);

        float Residu[N];
        for (int i = 0; i < N; i++) {
                Residu[i] = MatDD[i][i] * V4[i];
                for (int j = 0; j < N; j++) {
                        Residu[i] -= MatDD[i][j] * V4[j];
                }
        }
        float errorRelatiu = 0.0;
        for (int i = 0; i < N; i++){
                errorRelatiu += Residu[i] * Residu[i];
        }
        errorRelatiu = sqrt(errorRelatiu) / Magnitude(V3);
        printf("L'error relatiu és %f.\n",errorRelatiu);

        Jacobi(Mat, V3, V4, 1);
}
