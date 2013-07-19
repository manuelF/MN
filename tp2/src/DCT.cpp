#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <queue>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include <cassert>

#include <omp.h>
#include <math.h>

#include <random>

#define EPS 1e-8


#define GAUSSIAN_NOISE 0
#define ADDITIVE_NOISE 1
#define IMPULSE_NOISE 2

#define ZERO_FILTER 0
#define EXPONENTIAL_FILTER 1
#define MEDIAN_FILTER 3

#define IMPULSE_NOISE_DEFAULT 20
#define IMPULSE_NOISE_TICK 20


using namespace std;


const double max_m()
{
    return 255.0;
}

//Versiones para ecm de vector y de matriz
double ecm(const vector < double >&orig, const vector < double >&recup);
double ecm(const vector < vector < double >>&orig,
           const vector < vector < double >>&recup);


//ECM de Matriz
double ecm(const vector < vector < double > >&orig,
           const vector < vector < double > >&recup)
{
    double e = 0.0;

    if (orig.size() != recup.size())
    {
        cerr << "Error-ECM: distintas longitudes de matrices" << endl;
        exit(1);
    }
    unsigned int n = orig.size();

    for (unsigned int i = 0; i < n; i++)
    {
        double q = ecm(orig[i], recup[i]);
        e += (q * q);
    }
    e = e / (double) (n * n);
    return e;

}

//ECM de Vector
double ecm(const vector < double >&orig, const vector < double >&recup)
{
    double e = 0.0;

    if (orig.size() != recup.size())
    {
        cerr << "Error-ECM: distintas longitudes de vectores" << endl;
        exit(1);
    }
    unsigned int n = orig.size();

    for (unsigned int i = 0; i < n; i++)
    {
        double t = orig[i] - recup[i];
        e += (t * t);
    }
    e = e / (double) n;
    return e;
}


//PSNR de Matriz
double psnr(const vector < vector < double >>&orig,
            const vector < vector < double >>&recup)
{
    double _max = 255.0;
    return 10 * log10((_max * _max) / ecm(orig, recup));
}

//PSNR de Vector
double psnr(const vector < double >&orig, const vector < double >&recup)
{
    double _max = *max_element(orig.begin(), orig.end());

    return 10 * log10((_max * _max) / ecm(orig, recup));
}

vector < double >frecuencias;
vector < double >muestreo;
vector < double >lecturas;
vector < double >ruido;
vector < vector < double >>MSombrero;
vector < vector < double >>M;

//Dado B calculo MBM^t como pide en el apendice A1
vector < vector < double >>mbmt(const vector < vector < double > >&b)
{
    int n = M.size();
    //Tiene que ser cuadrada la matriz M, y compatible el tamaño de la matriz b
    assert((int) M[0].size() == n && (int) b.size() == n
           && (int) b[0].size() == n);
    // Creo temporal con todos ceros la matriz
    vector < vector < double >>aux(n, vector < double >(n, 0));

    //Multiplico las matrices
    #pragma omp parallel for
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            for (int t = 0; t < n; t++)
                aux[i][j] += M[i][t] * b[t][j];
        }
    }
    //Creo la matriz de resultado
    vector < vector < double >>sol(n, vector < double >(n, 0));
    //Aplico la multiplicacion para el otro lado con la transpuesta
    #pragma omp parallel for
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            for (int t = 0; t < n; t++)
            {
                sol[i][j] = aux[i][t] * M[j][t];	// (j,t) y no (t,j) porque es la transpuesta
            }
        }
    }
    //Devuelvo M x b x Mt
    return sol;
}

//Resuelvo el sistema de ecuaciones para deshacer la transformacion
vector < double >antitransformar(const vector < double >&y, const vector<vector<double>>& mat)
{
    int n = mat.size();
    //La matriz debe ser cuadrada
    assert((int) mat.size() == n && (int) mat[0].size() == n
           && (int) y.size() == n);
    // Me copio la matriz
    vector < vector < double >>sistema(mat);
    //Genero la matriz adjunta [mat|y]
    for (int i = 0; i < n; i++)
        sistema[i].push_back(y[i]);

    //Resuelvo el sistema
    for (int j = 0; j < n - 1; j++)
    {
        int mx = j;
        //Obtengo la minima fila
        for (int t = j + 1; t < n; t++)
            if (abs(sistema[mx][j]) < abs(sistema[t][j]))
                mx = t;
        //Permuto con la que estoy mirando (minimiza el error)
        swap(sistema[mx], sistema[j]);

        double m = sistema[j][j];

        //Hago un paso de Gauss
        for (int i = j + 1; i < n; i++)
        {
            double db = sistema[i][j] / m;
            for (int t = j; t < n + 1; t++)
                sistema[i][t] -= db * sistema[j][t];
        }
    }

    vector < double >x(n, 0);	// los valores de x antes de aplicar la permutacion
    for (int i = n - 1; i >= 0; i--)
    {
        double y_i = sistema[i][n];	// el que seria el numero en y en [mat|y]
        double a = 0;
        for (int j = i + 1; j < n; j++)
            a += sistema[i][j] * x[j];
        x[i] = (y_i - a) / sistema[i][i];
    }
    //Devuelvo el vector antitransformado
    return x;
}

//Aplicamos la formula del enunciado para generarlo
void generarMatrizDCT(int n)
{
    //Inicializamos en vacios los vectores globales
    frecuencias = vector < double >(n);
    muestreo = vector < double >(n);
    MSombrero = vector < vector < double >>(n, vector < double >(n));
    M = vector < vector < double >>(n, vector < double >(n));


    double k = (M_PI / (double) n);
    for (int i = 0; i < n; i++)
    {
        frecuencias[i] = i;
        muestreo[i] = k * (i + ((1. / 2.)));
    }
    double sq1n = sqrt(1.0 / n);
    double sq2n = sqrt(2.0 / n);
    double _max = max_m();
    double v = sq1n;
    //Generaro la matriz de la DCT
    #pragma omp parallel for
    for (int i = 0; i < n; i++)
    {
        double freq = frecuencias[i];
        for (int j = 0; j < n; j++)
        {
            if (i > 0)
                v = sq2n;
            double theta = v * cos(freq * muestreo[j]);
            MSombrero[i][j] = theta;
            M[i][j] = floor((_max * theta + 1) / 2);
        }
    }
}

//Transformar es una multiplicacion Matriz por Vector
vector < double >transformar(const vector < double >&l)
{
    int n = (int) l.size();
    vector < double >y(n);

    #pragma omp parallel for
    for (int i = 0; i < n; i++)
    {
        double q = 0.0;
        for (int j = 0; j < n; j++)
        {
            q += M[i][j] * l[j];
        }
        y[i] = q;

    }
    return y;
}

//Generamos un ruido sobre el vector pasado como parametro
void generarRuido(vector < double >&y, int imp)
{
    //Es un ruido blanco
    if (imp == GAUSSIAN_NOISE)
    {
        ruido = vector < double >(y.size());
        std::default_random_engine generator;
        std::normal_distribution < double >distribution(0.0, 10.0);
        for (int i = 0; i < (int) y.size(); i++)
        {
            double randval = distribution(generator);
            ruido[i] = randval;
            y[i] = y[i] + ruido[i];
        }
    }

    //Es un ruido impulsivo, muy angosto y fuerte
    if (imp == IMPULSE_NOISE)
    {
        int sign = 1;
        ruido = vector < double >(y.size());
        for (int i = 0; i < (int) y.size(); i++)
        {
            //Cada tanto (IMPULSE_NOISE_TICK), generamos un pico
            //(una vez para cada lado)
            if (i % IMPULSE_NOISE_TICK == 0)
            {
                if (i % 3 * IMPULSE_NOISE_TICK == 0)
                    sign *= -1;
                y[i] += sign * IMPULSE_NOISE_DEFAULT;
            }
        }
    }
    return;

}

//Aplicamos algun filtro a una tira de numeros
void filtrarRuido(vector < double >&y, int imp)
{
    const int n = (int) y.size();
    //Parametro: Contorno alrededor de un pico que vamos a apagar
    const int contorno = 20;
    //Parametro: Threshold de diferencia con el maximo para que se active el filtro
    const double diffmax = 2.0;
    //Parametro: Porcentaje de la muestra a partir de la cual vamos a mirar
    //(esto evita filtrar las componentes de muy baja frecuencia)
    const int onepercent = n / 100;
    const int startfilter = 5 * onepercent;

    //Vector de salida
    vector < double >replace(y);

    //Filtro que pone en ceros el contorno alrededor de los picos
    if (imp == ZERO_FILTER)
    {
        double mx = abs(y[startfilter]);
        for (int i = startfilter; i < n; i++)
            mx = max(mx, abs(y[i]));
        for (int i = startfilter; i < n; i++)
        {
            if (abs(y[i]) * diffmax > mx)
            {
                for (int j = max(startfilter, i - contorno);
                        j < min(n, i + contorno); j++)
                    replace[j] = 0;
            }
        }
    }
    //Filtro que disminuye exponencialmente alrededor de los picos
    //(preserva un poco mas la forma de la senal)
    if (imp == EXPONENTIAL_FILTER)
    {
        double mx = abs(y[startfilter]);
        for (int i = startfilter; i < n; i++)
            mx = max(mx, abs(y[i]));
        for (int i = startfilter; i < n; i++)
        {
            if (abs(y[i]) * diffmax > mx)
            {
                for (int j = max(0, i - contorno);
                        j < min(n, i + contorno); j++)
                {
                    double decay =
                        (1.0 /
                         pow(1.2, contorno - abs(j - i) + 1));
                    replace[j] = y[j] * decay;
                }

            }
        }
    }
    //Filtro que obtiene la mediana alrededor de los contornos
    if (imp == MEDIAN_FILTER)
    {
        vector < double >local_y;
        int qty = 4;
        for (int i = startfilter + qty; i < n - qty; i++)
        {

            local_y.clear();
            for (int j = i - qty; j <= i + qty; j++)
            {
                local_y.push_back(y[j]);
            }

            sort(local_y.begin(), local_y.end());
            replace[i] = local_y[qty];
        }

    }
    y=replace;
}

//Los filtros anteriores pero implementados para 2D
void filtrarRuido2D(vector < vector < double > >&y, int imp)
{
    const int n = (int) y.size();	// == y[0].size()
    //Asumimos imagen cuadrada
    assert((int)y[0].size()==n);

    const int contorno = 20;
    const double diffmax = 2.0;
    const int onepercent = n / 100;
    const int startfilter = 5 * onepercent;
    vector < vector < double >>replace(y);

    //Pone ceros alrededor de un pico
    if (imp == ZERO_FILTER)
    {
        double mx = 0.0;
        for (int i = startfilter; i < n; i++)
            for (int j = startfilter; j < n; j++)
                mx = max(mx, abs(y[i][j]));

        for (int i = startfilter; i < n; i++)
            for (int j = startfilter; j < n; j++)
            {
                if (abs(y[i][j]) * diffmax > mx)
                {
                    for (int ii = max(startfilter, i - contorno);
                            ii < min(n, i + contorno);
                            ii++)
                        for (int jj = max(startfilter, j - contorno);
                                jj < min(n,j + contorno); jj++)
                            replace[ii][jj] = 0;
                }
            }
    }
    //Hace decaimiento exponencial alrededor de un pico
    if (imp == EXPONENTIAL_FILTER)
    {
        double mx = 0.0;
        for (int i = startfilter; i < n; i++)
            for (int j = startfilter; j < n; j++)
                mx = max(mx, abs(y[i][j]));

        for (int i = startfilter; i < n; i++)
            for (int j = startfilter; j < n; j++)
            {
                if (abs(y[i][j]) * diffmax > mx)
                {
                    for (int ii = max(0, i - contorno);
                            ii <= min(n, i + contorno); ii++)
                        for (int jj = max(0, j - contorno);
                                jj <= min(n, j + contorno); jj++)
                        {
                            double decay = (1.0 / pow(1.2,
                                                      2 *
                                                      contorno -
                                                      abs(ii - i) -
                                                      abs(jj - j) +
                                                      1));
                            replace[ii][jj] = y[ii][jj] * decay;
                        }

                }
            }
    }
    //Reemplaza un punto por la mediana de todos los vecinos en un contorno
    if (imp == MEDIAN_FILTER)
    {
        vector < double >local_y;
        int qty = 4;
        for (int i = startfilter + qty; i < n - qty; i++)
            for (int j = startfilter + qty; j < n - qty; j++)
            {
                local_y.clear();
                for (int ii = i - qty; ii <= i + qty; ii++)
                    for (int jj = j - qty; jj <= j + qty; jj++)
                        local_y.push_back(y[ii][jj]);
                sort(local_y.begin(), local_y.end());
                replace[i][j] = local_y[(qty * qty) / 2];
            }
    }

    //Guardamos los cambios
    y=replace;

}


//Funcion que dumpea los contenidos de un vector
//en un formato compatible para que lo puedan tomar
//los scripts de gnuplot
void dump(string name, const vector < double >&src)
{
    int n = src.size();
    FILE *fileptr = fopen(name.c_str(), "w");
    for (int i = 0; i < n; i++)
        fprintf(fileptr, "%.6lf\n", src[i]);
    fclose(fileptr);
}

//Funcion cuasi principal para filtrar senales de una dimension
void procesar1D()
{
    //Tomamos n muestras, n la primer medicion
    int n;
    cin >> n;
    lecturas = vector < double >(n);
    double _max = 0.0;
    for (int i = 0; i < n; i++)
    {
        cin >> lecturas[i];
        _max = (_max > abs(lecturas[i]) ? _max : abs(lecturas[i]));
    }

    //Obtenemos la matriz de DCT apropiada para ese n

    generarMatrizDCT(n);
    //Duplicamos las lecturas para poder comparar el exito de la transformacion
    //despues
    vector < double >l(lecturas);

    //Generamos algun ruido sobre la senal

    //generarRuido(l,GAUSSIAN_NOISE);
    generarRuido(l, IMPULSE_NOISE);


    //Pasamos las senales al espacio DCT
    vector < double >q = transformar(lecturas);	// senal original

    vector < double >y = transformar(l);	//senal ruidosa

    dump("orig", q);
    dump("mod", y);

    //Aplicamos algun/os filtro/s sobre la senal ruidosa

    filtrarRuido(y, EXPONENTIAL_FILTER);
    //filtrarRuido(y,ZERO_FILTER );
    //filtrarRuido(y,MEDIAN_FILTER );

    dump("recovered", y);

    //Aplicamos la antitransformacion para obtener una reconstruccion
    //de la senal original
    vector < double >x = antitransformar(y,M);

    //Escribimos cual fue el PNSR calculado entre la senal original y la
    //reconstruida
    cerr << "PNSR: " << psnr(lecturas, x) << endl;
}

//Funcion auxiliar para transponer una matriz para obtener mejor performance
void transpose(vector < vector < double > >&mat)
{
    for (int i = 0; i < (int) mat.size(); i++)
    {
        for (int j = i; j < (int) mat[i].size(); j++)
        {
            swap(mat[i][j], mat[j][i]);
        }
    }
}
// A*ret=B
vector < vector < double >>resolverSistema(vector < vector < double > >&A, vector < vector < double > >&B)
{
    int n = A.size();		/// A[0].size() == B.size() == B[0].size()
    vector < vector < double >>ret(n);
    transpose(B);
    for (int i = 0; i < n; i++)
        ret[i] = antitransformar(B[i],A);
    transpose(B);
    transpose(ret);
    return ret;
}

vector < vector < double >>antitransformar2D(vector < vector < double > >&y)	/// mat * ret * mat^t= y
{
    /**
     * Primero hacemos mat*x = y, y despues ret*mat^t = x que es o mismo que mat*ret^t = x^t
     */
    vector < vector < double >>mat = M;
    vector < vector < double >>x = resolverSistema(mat, y);
    transpose(x);
    vector < vector < double >>ret = resolverSistema(mat, x);
    transpose(ret);
    return ret;
}

/**
 *
 * tranformar2D es M mat M^t. Primero hago mat M^t y despues hago M por el resultado de eso
 **/

vector < vector < double >>transformar2D(const vector < vector <
        double > >&mat)
{
    int n = mat.size();
    vector < vector < double >>y (n, vector < double >(n, 0)),
           res(n, vector < double>(n, 0));

    #pragma omp parallel for
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
        {
            for (int t = 0; t < n; t++)
                y[i][j] += mat[i][t] * M[j][t];	///M[j][t] porque es M^t
        }
    #pragma omp parallel for
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
        {
            for (int t = 0; t < n; t++)
                res[i][j] += M[i][t] * y[t][j];
        }
    return res;
}

//Funcion cuasi principal que ordena los filtros a las imagenes
void procesar2D()
{
    //Tiene que matchear el magic token
    string magic;
    cin >> magic;

    if (magic != "P5")
    {
        cerr << "ERROR: PGM tipo " << magic <<
        " no implementado, solo P5" << endl;
        exit(1);
    }
    //Levantamos la meta data
    int width, height;
    int grayscale;
    int scanf_res = scanf("%d %d\n%d\n", &width, &height, &grayscale);
    ///recordemos que width == height siempre
    if (scanf_res != 3)
    {
        cerr << "ERROR: Leida no concuerda con esperado" << endl;
    }

    vector < vector < double >>_img(height, vector < double >(width));
    //Levantamos toda la imagen
    vector < double >todo(width * height);
    for (int j = 0; j < height; j++)
    {
        for (int i = 0; i < width; i++)
        {
            unsigned char read;
            if (scanf("%c", &read) != 1)
            {
                cerr << "ERROR: PGM no se puede leer correctamente" << endl;
            }
            _img[j][i] = (double) read;
            todo[j * width + i] = (double) read;
        }
    }
    //Generamos la matriz de transformacion
    generarMatrizDCT(width);
    //Abrimos el archivo de salida para guarda la imagen transformada
    FILE *f = fopen("imgMod.pgm", "w");
    fprintf(f, "P5\n%d %d\n%d\n", width, height, grayscale);

    //Guardamos la original para obtener el PSNR
    vector < vector < double >>vec = _img;

    //Agregamos un ruido fila a fila
    for (int j = 0; j < height; j++)
        generarRuido(vec[j], IMPULSE_NOISE);

    //Transformamos la senal ruidosa al espacio DCT
    vector < vector < double >>vec2 = transformar2D(vec);	//transformada

    //Aplicamos un filtro 2D usando bloques
    filtrarRuido2D(vec2, EXPONENTIAL_FILTER);

    //Resolvemos el sistema para 2D
    vector < vector < double >>out = antitransformar2D(vec2);

    //Obtenemos el PSNR de lo bien que salio esta pasada
    double e = psnr(_img, out);
    cerr << "PSNR: " << e << endl;
    
    //Guardamos la imagen
    for (int j = 0; j < height; j++)
    {
        for (int i = 0; i < width; i++)
        {
            fprintf(f, "%c", (unsigned char) out[j][i]);
        }
    }
    //Cerramos para flushear
    fclose(f);

    //Fin de proceso
    return;

}
//Imprime como se usa este ejecutable
void usage()
{
    cerr << "Uso: ./DCT <imp>" << endl;
    cerr << "donde imp = (opcional) imp = 0 | Lo que se lee por stdin es una imagen  "<< endl;
    cerr << "                       imp = 1 | Lo que se lee por stdin es una senal 1D"<< endl;

    return;
}
//Funcion principal, delega el filtrado correcto
int main(int argc, char *argv[])
{
    switch (argc)
    {
        //Sin parametros, asumimos que leemos por stdin una senal 1D
    case 1:
        procesar1D();
        break;
    case 2:
        //Si es cero, asumimos que es una senal por stdin
        if (atoi(argv[1]) == 0)
            procesar1D();
        //Con un parametro, si es 1 es una imagen por stdin,
        if (atoi(argv[1]) == 1)
            procesar2D();
        //Otro valor, usage y out
        if (atoi(argv[1]) > 1)
            usage();
        break;
    default:
        //Otra cantidad, mandamos usage
        usage();
        break;

    }
    return 0;

}
