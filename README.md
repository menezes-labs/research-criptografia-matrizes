# matrixcalculus1000

[![CI](https://github.com/menezesx2k26-byte/matrixcalculus1000/actions/workflows/ci.yml/badge.svg)](https://github.com/menezesx2k26-byte/matrixcalculus1000/actions/workflows/ci.yml)

Projeto em **Rust** para gerar, analisar e medir o custo de matrizes reais densas quadradas de até **1000×1000**, com foco em detectar matrizes **sem inversa** sem desperdiçar trabalho calculando a inversa completa.

A decisão matemática central é:

```text
A é invertível  <=>  rank(A) = n
A é singular    <=>  rank(A) < n
```

O posto é obtido por **eliminação de Gauss com pivoteamento parcial**. A implementação usa armazenamento contíguo (`Vec<f64>`), tolerância numérica relativa à escala da matriz e paralelização com Rayon nas eliminações de matrizes maiores.

## O que o projeto faz

- aceita matrizes quadradas de ordem `1..=1000`;
- detecta singularidade numérica e calcula posto e nulidade;
- não modifica a matriz original durante a análise;
- usa pivoteamento parcial para melhorar estabilidade numérica;
- calcula automaticamente uma tolerância baseada em `epsilon`, ordem e norma infinito;
- permite sobrescrever a tolerância pela CLI;
- paraleliza a eliminação das linhas abaixo do pivô em matrizes grandes;
- armazena os elementos em um único bloco contínuo de memória;
- gera casos singulares deterministicamente;
- gera casos garantidamente invertíveis por **dominância diagonal estrita**;
- gera matrizes aleatórias para experimentos;
- oferece repetição de execuções e resumo de benchmark;
- exporta resultados em CSV;
- possui testes automatizados e CI no GitHub Actions.

## Por que não calcular a inversa?

Se a pergunta é apenas “esta matriz possui inversa?”, calcular `A⁻¹` inteiro é trabalho desnecessário.

Para `A ∈ Mₙ(ℝ)`, as seguintes afirmações são equivalentes:

```text
det(A) ≠ 0
rank(A) = n
ker(A) = {0}
A é invertível
```

A implementação acompanha os pivôs durante a eliminação. Se não existem `n` pivôs numericamente válidos, a matriz possui posto deficiente e é classificada como singular.

## Estrutura

```text
.
├── Cargo.toml
├── rust-toolchain.toml
├── src
│   ├── lib.rs
│   ├── matrix.rs       # estrutura densa + eliminação + análise
│   ├── generator.rs    # casos singular, invertível e aleatório
│   └── main.rs         # CLI + benchmark
├── tests
│   └── matrix_tests.rs
└── .github
    └── workflows
        └── ci.yml
```

## Requisitos

- Rust stable
- Cargo

O arquivo `rust-toolchain.toml` mantém o projeto na toolchain stable e inclui `rustfmt` e `clippy`.

## Executar

### Matriz singular 1000×1000

```bash
cargo run --release -- --size 1000 --kind singular
```

### Matriz garantidamente invertível 1000×1000

```bash
cargo run --release -- --size 1000 --kind invertible
```

### Repetir cinco vezes para medir desempenho

```bash
cargo run --release -- --size 1000 --kind singular --repeat 5
```

### Fixar quantidade de threads

```bash
cargo run --release -- --size 1000 --kind singular --threads 8
```

### Exportar benchmark como CSV

```bash
cargo run --release -- --size 500 --kind singular --repeat 10 --csv > benchmark.csv
```

### Imprimir matriz pequena

```bash
cargo run --release -- --size 5 --kind singular --print
```

A impressão é limitada a matrizes de até `20×20` para evitar despejar milhões de números no terminal.

## Opções da CLI

| Opção | Função | Padrão |
|---|---|---:|
| `--size N` | ordem da matriz, `1..=1000` | `5` |
| `--kind KIND` | `singular`, `invertible` ou `random` | `singular` |
| `--seed N` | seed determinística `u64` | `42` |
| `--repeat N` | número de execuções | `1` |
| `--tolerance VALUE` | tolerância numérica ou `auto` | `auto` |
| `--threads N` | workers Rayon | CPUs lógicas |
| `--print` | imprime matrizes até `20×20` | desligado |
| `--csv` | saída tabular para benchmark | desligado |

## Geração dos casos de teste

### `singular`

Para `n ≥ 3`, a última linha é construída como combinação linear das duas primeiras:

```text
Lₙ = 1.25 L₁ - 0.75 L₂
```

Logo as linhas são linearmente dependentes e a matriz é singular por construção.

Para `n = 2`, a segunda linha é o dobro da primeira. Para `n = 1`, o único elemento é zero.

### `invertible`

Cada linha é construída de forma que:

```text
|aᵢᵢ| > Σ |aᵢⱼ|,  j ≠ i
```

Portanto a matriz é estritamente diagonalmente dominante por linhas e, pelo teorema de Levy–Desplanques, é não singular.

### `random`

Todos os elementos são gerados no intervalo `[-1, 1)`. Nesse modo não existe expectativa prévia sobre a classificação.

A geração usa `SplitMix64` interno para manter os experimentos determinísticos sem adicionar outra dependência externa.

## Algoritmo

A análise segue, em alto nível:

1. copia a matriz para preservar a entrada;
2. calcula a tolerância automática, se nenhuma for fornecida;
3. percorre as colunas procurando o maior pivô em módulo;
4. troca linhas quando necessário;
5. elimina os elementos abaixo do pivô;
6. paraleliza as linhas independentes em matrizes suficientemente grandes;
7. conta os pivôs encontrados;
8. retorna posto, nulidade, singularidade e estatísticas auxiliares.

A análise também acumula o sinal do determinante e `ln|det(A)|` quando a matriz é classificada como invertível. O log é usado porque multiplicar centenas de pivôs diretamente pode causar overflow ou underflow em `f64`.

## Complexidade

Para uma matriz densa `n×n`:

```text
tempo:   O(n³)
memória: O(n²)
```

Uma matriz `1000×1000` contém:

```text
1.000.000 elementos
```

Em `f64`, apenas os dados ocupam:

```text
1.000.000 × 8 bytes = 8.000.000 bytes ≈ 7,63 MiB
```

Durante a análise existe uma cópia de trabalho da matriz, então o consumo principal fica perto de duas matrizes densas, além de buffers pequenos e overhead do runtime.

## Tolerância e ponto flutuante

Este projeto trabalha sobre aproximações em `f64`. Portanto “singular” significa **numericamente singular segundo a tolerância escolhida**.

A tolerância automática usa a ordem da matriz, a norma infinito e `f64::EPSILON`. Isso evita depender de um valor absoluto fixo que funcionaria mal quando a escala dos elementos muda muito.

Se o objetivo for álgebra exata, inteiros, racionais ou aritmética modular, o domínio matemático deve ser representado explicitamente em vez de usar `f64`.

## Atenção para criptografia modular

Uma matriz pode ser invertível sobre `ℝ` e não ser invertível módulo `m`.

Em `ℤₘ`, o critério é:

```text
gcd(det(A), m) = 1
```

Por exemplo, com módulo `26`, uma matriz cujo determinante seja `2` é invertível sobre os reais, mas não em `ℤ₂₆`, porque `gcd(2, 26) = 2`.

Por isso este núcleo real não deve ser usado como substituto direto de um teste de invertibilidade para Hill Cipher ou outra construção modular.

## Testes

```bash
cargo test
```

A suíte cobre:

- matriz singular conhecida;
- identidade invertível;
- geração singular;
- geração estritamente diagonalmente dominante;
- alocação estrutural de `1000×1000`;
- rejeição de ordem acima do limite;
- rejeição de `NaN`/infinito;
- rejeição de tolerância inválida.

## Qualidade e CI

O GitHub Actions executa em cada push/PR:

```bash
cargo fmt --all -- --check
cargo clippy --all-targets --all-features
cargo test --all
cargo build --release
```

Nenhum número de benchmark é publicado no README sem ter sido realmente medido. Use `--repeat` e `--csv` para produzir resultados reproduzíveis na máquina desejada.

## Licença

MIT.
