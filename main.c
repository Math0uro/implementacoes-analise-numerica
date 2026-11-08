#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_ITERACOES 100
#define ER_RAIZES_NEWTON 0.000001
#define ER_OUTROS 0.00001
#define N 3

typedef struct {
    int convergiu;
    int iteracoes;
    double erro;
    double valor;
} Resultado;

static double erro_relativo(double novo, double anterior) {
    double distancia = fabs(novo - anterior);
    return novo != 0.0 ? distancia / fabs(novo) : distancia;
}

static void mostrar_resultado_raiz(Resultado r, int erro_absoluto) {
    printf("Situacao: %s\n", r.convergiu ? "Convergiu" : "Nao convergiu");
    printf("Iteracoes: %d / %d\n", r.iteracoes, MAX_ITERACOES);
    if (isfinite(r.valor)) printf("Raiz aproximada: %.10f\n", r.valor);
    if (isfinite(r.erro)) printf("Erro %s: %.9g\n", erro_absoluto ? "absoluto estimado" : "relativo", r.erro);
}

static double funcao_newton(double x) { return x - 4.0 + log(x); }
static double derivada_newton(double x) { return 1.0 + 1.0 / x; }
static double funcao_ponto_fixo(double x) { return -x*x*x/3.0 + 3.0*x - 1.0; }
static double iteracao_ponto_fixo(double x) { return x*x*x/9.0 + 1.0/3.0; }
static double funcao_bissecao(double x) { return exp(x) - x*x - 3.0*x; }

static Resultado newton(void) {
    Resultado r = {0, 0, INFINITY, 2.0};
    for (int k = 1; k <= MAX_ITERACOES; ++k) {
        double x = r.valor;
        if (x <= 0 || !isfinite(x)) break;
        double derivada = derivada_newton(x);
        if (!isfinite(derivada) || fabs(derivada) < 1e-14) break;
        double proximo = x - funcao_newton(x) / derivada;
        if (!isfinite(proximo) || proximo <= 0) break;
        r.iteracoes = k;
        r.erro = erro_relativo(proximo, x);
        r.valor = proximo;
        if (r.erro < ER_RAIZES_NEWTON) { r.convergiu = 1; break; }
    }
    return r;
}

static Resultado ponto_fixo(void) {
    Resultado r = {0, 0, INFINITY, 0.0};
    for (int k = 1; k <= MAX_ITERACOES; ++k) {
        double proximo = iteracao_ponto_fixo(r.valor);
        if (!isfinite(proximo) || proximo < 0 || proximo > 1) break;
        double erro = erro_relativo(proximo, r.valor);
        r.iteracoes = k;
        r.erro = erro;
        r.valor = proximo;
        if (erro < ER_OUTROS) { r.convergiu = 1; break; }
    }
    return r;
}

static Resultado bissecao(void) {
    double a = -4.0, b = -3.0, fa = funcao_bissecao(a);
    Resultado r = {0, 0, INFINITY, NAN};
    if (fa * funcao_bissecao(b) >= 0) return r;
    for (int k = 1; k <= MAX_ITERACOES; ++k) {
        double meio = (a + b) / 2.0;
        double fm = funcao_bissecao(meio);
        r.iteracoes = k;
        r.valor = meio;
        r.erro = (b - a) / 2.0;
        if (fm == 0 || (b - a) < ER_OUTROS) { r.convergiu = 1; break; }
        if (fa * fm < 0) b = meio;
        else { a = meio; fa = fm; }
    }
    return r;
}

static double norma_maxima(const double v[N]) {
    double maior = 0.0;
    for (int i = 0; i < N; ++i)
        if (fabs(v[i]) > maior) maior = fabs(v[i]);
    return maior;
}

static void verificar_criterios(const double A[N][N]) {
    double beta[N] = {0};
    int linhas = 1;
    int sassenfeld = 1;
    for (int i = 0; i < N; ++i) {
        double soma = 0.0, soma_sassenfeld = 0.0;
        for (int j = 0; j < N; ++j) {
            if (i == j) continue;
            soma += fabs(A[i][j]);
            soma_sassenfeld += fabs(A[i][j]) * (j < i ? beta[j] : 1.0);
        }
        if (fabs(A[i][i]) <= soma) linhas = 0;
        beta[i] = soma_sassenfeld / fabs(A[i][i]);
        if (beta[i] >= 1.0) sassenfeld = 0;
    }
    printf("Criterio das linhas: %s\n", linhas ? "Satisfeito" : "Nao satisfeito");
    printf("Criterio de Sassenfeld: %s\n", sassenfeld ? "Satisfeito" : "Nao satisfeito");
}

static Resultado gauss_seidel(double x[N]) {
    /* Mesmas equacoes do enunciado, na ordem 3, 2, 1. */
    const double A[N][N] = {{2, 2, 1}, {1, 1, 1}, {1, 2, -2}};
    const double b[N] = {5, 2, 7};
    Resultado r = {0, 0, INFINITY, NAN};

    puts("Sistema original: (1, 2, -2), (1, 1, 1), (2, 2, 1)");
    puts("Ordem original:");
    const double original[N][N] = {{1, 2, -2}, {1, 1, 1}, {2, 2, 1}};
    verificar_criterios(original);
    puts("Equacoes reordenadas: 3, 2, 1");
    verificar_criterios(A);
    puts("Os criterios nao garantem convergencia; vamos verificar pelas iteracoes.");

    for (int k = 1; k <= MAX_ITERACOES; ++k) {
        double anterior[N];
        for (int i = 0; i < N; ++i) anterior[i] = x[i];
        for (int i = 0; i < N; ++i) {
            double soma = 0.0;
            for (int j = 0; j < N; ++j)
                if (i != j) soma += A[i][j] * x[j];
            x[i] = (b[i] - soma) / A[i][i];
        }
        double diferenca[N];
        for (int i = 0; i < N; ++i) diferenca[i] = x[i] - anterior[i];
        double denominador = norma_maxima(x);
        r.erro = denominador > 0 ? norma_maxima(diferenca) / denominador
                                : norma_maxima(diferenca);
        r.iteracoes = k;
        if (!isfinite(r.erro)) break;
        if (r.erro < ER_OUTROS) { r.convergiu = 1; break; }
    }
    return r;
}

static int executar(const char *metodo) {
    if (strcmp(metodo, "newton") == 0) {
        puts("NEWTON | f(x) = x - 4 + ln(x)");
        puts("x0 = 2 | Er < 0.000001");
        Resultado r = newton();
        mostrar_resultado_raiz(r, 0);
        return !r.convergiu;
    }
    if (strcmp(metodo, "ponto-fixo") == 0) {
        puts("PONTO FIXO | f(x) = -x^3/3 + 3x - 1");
        puts("phi(x) = x^3/9 + 1/3 | I = [0,1] | x0 = 0 | Er < 0.00001");
        puts("Convergencia: |phi'(x)| < 1 em [0,1], em todo [0,1].");
        Resultado r = ponto_fixo();
        mostrar_resultado_raiz(r, 0);
        if (isfinite(r.valor)) printf("f(x) = %.9g\n", funcao_ponto_fixo(r.valor));
        return !r.convergiu;
    }
    if (strcmp(metodo, "bissecao") == 0) {
        puts("BISSECAO | f(x) = e^x - x^2 - 3x");
        puts("I = [-4,-3] | Precisao < 0.00001");
        Resultado r = bissecao();
        printf("Comprimento do intervalo: %.9g\n", 2*r.erro);
        mostrar_resultado_raiz(r, 1);
        return !r.convergiu;
    }
    if (strcmp(metodo, "gauss-seidel") == 0) {
        puts("GAUSS-SEIDEL | X0 = (0,0,0) | Er < 0.00001");
        double x[N] = {0.0, 0.0, 0.0};
        Resultado r = gauss_seidel(x);
        printf("Situacao: %s\n", r.convergiu ? "Convergiu" : "Nao convergiu");
        printf("Iteracoes: %d / %d\n", r.iteracoes, MAX_ITERACOES);
        printf("x1 = %.8f\nx2 = %.8f\nx3 = %.8f\n", x[0], x[1], x[2]);
        if (isfinite(r.erro)) printf("Erro relativo: %.9g\n", r.erro);
        return !r.convergiu;
    }
    fprintf(stderr, "Metodo desconhecido: %s\n", metodo);
    return 2;
}

int main(int argc, char **argv) {
    if (argc == 2) return executar(argv[1]);
    if (argc != 1) {
        fprintf(stderr, "Uso: %s [newton|ponto-fixo|bissecao|gauss-seidel]\n", argv[0]);
        return 2;
    }
    const char *metodos[] = {"newton", "ponto-fixo", "bissecao", "gauss-seidel"};
    int falhas = 0;
    for (int i = 0; i < 4; ++i) {
        if (i) puts("\n------------------------------");
        falhas |= executar(metodos[i]);
    }
    return falhas ? EXIT_FAILURE : EXIT_SUCCESS;
}
