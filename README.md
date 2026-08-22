# matrixcalculus1000

[![CI](https://github.com/menezesx2k26-byte/matrixcalculus1000/actions/workflows/ci.yml/badge.svg)](https://github.com/menezesx2k26-byte/matrixcalculus1000/actions/workflows/ci.yml)

Projeto em **Rust** para gerar, analisar e medir o custo de matrizes reais densas quadradas de até **1000×1000**, com foco em detectar matrizes **sem inversa** sem calcular a inversa completa.

A ideia matemática central é:

```text
A é invertível  <=>  rank(A) = n
A é singular    <=>  rank(A) < n
```

O posto é obtido por **eliminação de Gauss com pivoteamento parcial**. A implementação usa armazenamento contíguo (`Vec<f64>`), tolerância numérica relativa à escala da matriz e paralelização com Rayon nas eliminações de matrizes maiores.

## Recursos

- matrizes quadradas de ordem `1..=1000`;
- posto, nulidade e detecção de singularidade;
- matriz original preservada durante a análise;
- pivoteamento parcial;
- tolerância automática ou configurável;
- eliminação paralela com Rayon para matrizes grandes;
- memória contígua em vez de `Vec<Vec<f64>>`;
- gerador determinístico de matrizes singulares;
- gerador de matrizes garantidamente invertíveis por dominância diagonal estrita;
- gerador aleatório para experimentos;
- benchmark por repetição;
- saída CSV;
- testes automatizados;
- CI com `cargo check`, Clippy, testes e build release.

## Por que não calcular a inversa?

Se a pergunta é apenas “esta matriz possui inversa?”, calcular `A⁻¹` inteiro é trabalho desnecessário.

Para `A ∈ Mₙ(ℝ)`:

```text
det(A) ≠ 0
rank(A) = n
ker(A) = {0}
A é invertível
```

são condições equivalentes. O programa conta pivôs numericamente válidos durante a eliminação. Se não existem `n` pivôs, a matriz é classificada como singular.

## Estrutura

```text
.
├── Cargo.toml
├── rust-toolchain.toml
├── src
│   ├── lib.rs
│   ├── matrix.rs       # matriz densa + eliminação + análise
│   ├── generator.rs    # singular, invertível e aleatória
│   └── main.rs         # CLI + benchmark
├── tests
│   └── matrix_tests.rs
└── .github
    └── workflows
        └── ci.yml
```

## Executar

Requer Rust stable e Cargo.

### Singular 1000×1000

```bash
cargo run --release -- --size 1000 --kind singular
```

### Garantidamente invertível 1000×1000

```bash
cargo run --release -- --size 1000 --kind invertible
```

### Benchmark com cinco execuções

```bash
cargo run --release -- --size 1000 --kind singular --repeat 5
```

### Fixar quantidade de threads

```bash
cargo run --release -- --size 1000 --kind singular --threads 8
```

### Exportar CSV

```bash
cargo run --release -- --size 500 --kind singular --repeat 10 --csv > benchmark.csv
```

### Mostrar uma matriz pequena

```bash
cargo run --release -- --size 5 --kind singular --print
```

A impressão é limitada a `20×20`.

## CLI

| Opção | Função | Padrão |
|---|---|---:|
| `--size N` | ordem, `1..=1000` | `5` |
| `--kind KIND` | `singular`, `invertible`, `random` | `singular` |
| `--seed N` | seed determinística `u64` | `42` |
| `--repeat N` | repetições | `1` |
| `--tolerance VALUE` | tolerância ou `auto` | `auto` |
| `--threads N` | workers Rayon | CPUs lógicas |
| `--print` | imprime matriz até `20×20` | desligado |
| `--csv` | saída para benchmark | desligado |

## Como os casos são construídos

### `singular`

Para `n ≥ 3`, a última linha é uma combinação linear das duas primeiras:

```text
Lₙ = 1.25 L₁ - 0.75 L₂
```

Logo as linhas são linearmente dependentes. Para `n = 2`, `L₂ = 2L₁`; para `n = 1`, o único elemento é zero.

### `invertible`

Cada linha satisfaz:

```text
|aᵢᵢ| > Σ |aᵢⱼ|,  j ≠ i
```

A matriz é estritamente diagonalmente dominante por linhas e, pelo teorema de Levy–Desplanques, é não singular.

### `random`

Elementos pseudoaleatórios em `[-1, 1)`, usando `SplitMix64` interno. O mesmo seed reproduz a mesma matriz.

## Algoritmo

1. valida ordem e valores;
2. calcula a norma infinito;
3. define tolerância automática se necessário;
4. copia a matriz;
5. escolhe o maior pivô em módulo na coluna;
6. troca linhas;
7. elimina as linhas abaixo do pivô;
8. paraleliza essas linhas quando o problema é grande;
9. conta pivôs;
10. retorna posto, nulidade e singularidade.

Para matrizes invertíveis, a análise também acumula o sinal do determinante e `ln|det(A)|`. O log evita overflow/underflow causado pelo produto direto de centenas de pivôs.

## Complexidade

Para matriz densa `n×n`:

```text
tempo:   O(n³)
memória: O(n²)
```

Uma `1000×1000` possui 1.000.000 de elementos. Em `f64`, os dados ocupam 8.000.000 bytes, aproximadamente **7,63 MiB**. Durante a análise existe uma cópia de trabalho, então o consumo principal fica próximo de duas matrizes densas, além de um pequeno buffer de pivô e do runtime.

## Precisão numérica

O domínio deste núcleo é `f64`. Assim, “singular” significa **numericamente singular segundo a tolerância escolhida**.

A tolerância automática depende de:

```text
f64::EPSILON × ordem × norma_infinito × fator_de_segurança
```

Isso é mais robusto que um `EPS = 1e-10` fixo para matrizes em escalas muito diferentes.

Para álgebra exata, racionais, inteiros ou aritmética modular, o domínio deve ser modelado explicitamente.

## Criptografia modular

Uma matriz pode ser invertível sobre `ℝ` e não ser invertível módulo `m`.

Em `ℤₘ`, o critério é:

```text
gcd(det(A), m) = 1
```

Exemplo: `det(A) = 2` é diferente de zero nos reais, mas não é unidade em `ℤ₂₆`, pois `gcd(2, 26) = 2`.

Portanto este núcleo real não deve ser usado como substituto direto de um teste de invertibilidade para Hill Cipher ou outra construção modular.

## Testes

```bash
cargo test
```

A suíte cobre matriz singular conhecida, identidade, geração singular, geração estritamente diagonalmente dominante, alocação estrutural `1000×1000`, ordem acima do limite, `NaN`/infinito e tolerância inválida.

## CI

A pipeline do GitHub Actions executa em push e pull request:

```bash
cargo check --all-targets --all-features
cargo clippy --all-targets --all-features
cargo test --all
cargo build --release
```

O README não publica números de benchmark inventados. Para medir uma máquina real, use `--repeat` e `--csv`.

## Licença

MIT.
