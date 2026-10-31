# Solver CFD 2D em MATLAB

Solver de volumes finitos para escoamento incompressível, laminar e em regime permanente, em malha deslocada (*staggered*) uniforme. O projeto compara os algoritmos de acoplamento pressão-velocidade SIMPLE, SIMPLEC e PISO em quatro escoamentos e com sete esquemas convectivos.

Não há dependência de toolbox. O código foi executado no MATLAB R2026a.

## Como rodar

Tudo passa por `cfd_run.m`. Execute a partir da raiz do repositório:

```bash
# Simulação
matlab -batch "cfd_run LidDrivenCavity algoritmo SIMPLEC esquema quick Re 400 Nx 64 Ny 64"

# Valores padrão do escoamento
matlab -batch "cfd_run FlowOverStep"

# Suíte comparativa
matlab -batch "cfd_run compare estudos 135 fast 1"

# Testes
matlab -batch "cfd_run test"
```

Dentro do MATLAB, digite o mesmo comando (`cfd_run FlowOverStep algoritmo PISO`) ou chame como função (`cfd_run("FlowOverStep", "algoritmo", "PISO", "Re", 50)`). O nome do escoamento pode ser abreviado por um prefixo, como `Lid`. `cfd_run` sem argumentos mostra a ajuda.

### Opções da simulação

| Opção | Valores | Padrão |
|---|---|---|
| `algoritmo` | `SIMPLE`, `SIMPLEC`, `PISO` | `PISO` |
| `esquema` | `upwind1`, `upwind2`, `central`, `powerlaw`, `quick`, `muscl`, `hric` | `upwind1` |

Também é possível sobrescrever qualquer parâmetro do escoamento. Um nome que não existe naquele escoamento gera erro com a lista dos nomes válidos.

| Parâmetro | Couette-Poiseuille | Cavidade com tampa | Cavidade aberta | Degrau |
|---|---|---|---|---|
| Nome do escoamento | `CouettePoiseuille` | `LidDrivenCavity` | `FlowOverCavity` | `FlowOverStep` |
| Geometria [m] | `Lx`=5, `Ly`=1 | `L`=1 | `Ly_canal`=1, `L_cav`=1, `D_cav`=0.5, `Lx_ante`=1, `Lx_post`=1 | `H_step`=0.5, `Ly_canal`=1, `Lx_ante`=1, `Lx_post`=4 |
| Viscosidade | `mu`=0.1 | `Re`=1000 | `Re`=147 (em `Ly_canal`) | `Re`=100 (em `H_step`) |
| `U` [m/s] | 2 (placa superior) | 1 (tampa) | 1 (média na entrada) | 1 (média na entrada) |
| `rho` [kg/m³] | 1 | 1 | 1 | 1 |
| Gradiente | `G`=4 | não se aplica | não se aplica | não se aplica |
| `Nx` × `Ny` | 50 × 20 | 64 × 64 | 60 × 30 | 80 × 30 |
| `alpha_u`, `alpha_p` | 0.8, 0.8 | 0.5, 0.5 | 0.5, 0.5 | 0.5, 0.5 |
| `tol` | 1e-8 | 1e-6 | 1e-6 | 1e-6 |
| `max_iter` | 5000 | 5000 | 8000 | 8000 |

### Opções de `compare`

| Opção | Significado | Padrão |
|---|---|---|
| `estudos` | Um dígito por estudo, entre 1, 3, 4 e 5 | `135` |
| `fast` | `1` usa malhas menores e tolerâncias mais frouxas | `0` |

| Estudo | Conteúdo |
|---|---|
| 1 | Precisão: esquema × algoritmo, com erros L2/Linf e tempo de CPU |
| 3 | Sub-relaxação: mapa `alpha_u` × `alpha_p` × algoritmo |
| 4 | Sensibilidade do Couette a `mu`, `G` e `U` |
| 5 | Robustez: os 4 escoamentos × algoritmo |

Não existe estudo 2. O estudo 4 estava desativado no script original e por isso fica fora do padrão.

## Saídas

Cada simulação cria `results/<escoamento>_<algoritmo>_<esquema>_<data_hora>/`, sem sobrescrever execuções anteriores, com:

- `convergencia.png`, `linhas_corrente.png`, `campo_u.png`, `campo_v.png`, `magnitude_vetores.png` e `pressao.png`;
- `perfil_u.png` no Couette, e `perfil_u.png` com `perfil_v.png` na cavidade com tampa;
- `resultado.mat`, com o caso `c`, o resultado `r`, `algoritmo` e `esquema`.

A suíte grava em `results/comparacao_parte_1/{mat,figuras}`, e aqui uma nova execução sobrescreve os arquivos anteriores. A pasta `results/` não é versionada.

## Estrutura

```
cfd_run.m            entrada única
src/
  cases/             build_case: geometria, malha, parâmetros e campos iniciais
  discretization/    coeficientes, esquemas convectivos, fluxos e contornos
  solver/            laço SIMPLE/SIMPLEC/PISO, preditor, correção de pressão, resíduo
  postprocess/       interpolação para os centros, erro, figuras
  studies/           suíte comparativa
tests/               verificações do solver
```

O fluxo principal é `build_case` → `run_case` → `plot_results`. Em cada iteração, `run_case` chama:

1. `momentum_predictor`, que resolve `u*` e `v*` com `mass_flux`, `espac_disc`, `boundary_conditions` e `pressure_interp`;
2. `mass_residual`, para o resíduo de massa;
3. `compute_d_field` e `pressure_correc_field`, para obter `p'`;
4. a correção de velocidade e pressão. O PISO acrescenta aqui o 2º corretor;
5. `convergence`, que avalia a norma do resíduo do preditor.

## Método numérico

- **Malha:** pressão nos centros dos VCs; `u` e `v` nas faces. `u` tem dimensão `(Nx-1, Ny)` e `v`, `(Nx, Ny-1)`.
- **Sólidos:** a cavidade aberta e o degrau usam a malha do domínio retangular completo, com máscara `is_wall` nos VCs sólidos.
- **SIMPLE:** `d = A/(aP/alpha_u)`.
- **SIMPLEC:** `d = A/(aP/alpha_u - Σanb)`, voltando à fórmula do SIMPLE no nó em que o denominador fica abaixo de 10% de `aP/alpha_u`.
- **PISO:** um preditor e dois corretores por iteração externa, com a correção de vizinhança no segundo.
- **Equação de `p'`:** Gauss-Seidel com até 200 varreduras. Só a cavidade com tampa precisa fixar uma referência de pressão.
- **Esquemas:** upwind de 1ª ordem, diferenças centrais e power law entram direto nos coeficientes. Upwind de 2ª ordem (LUD), QUICK, MUSCL (limitador de van Albada) e HRIC usam correção diferida sobre a matriz de upwind.

### Condições de contorno

| Escoamento | Entrada | Saída | Paredes |
|---|---|---|---|
| Couette-Poiseuille | Perfil analítico | Gradiente nulo | Inferior fixa, superior a `U` |
| Cavidade com tampa | Não há | Não há | Tampa a `U`, demais fixas |
| Cavidade aberta | Poiseuille no canal | Gradiente nulo | Fixas |
| Degrau | Poiseuille no canal | Gradiente nulo | Fixas |

No Couette-Poiseuille, `G` entra apenas pelo perfil de entrada e pela solução analítica. O gradiente não é um termo fonte da equação de momento.

## Testes

`cfd_run test` executa `tests/test_solver.m`, que verifica:

- a convergência dos três algoritmos no Couette e a ordem do erro ao refinar em y;
- a conservação de massa e a recirculação na cavidade aberta com QUICK;
- a rejeição de parâmetro e de algoritmo inexistentes.
