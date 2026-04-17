---
title: "F-тест (критерий Фишера)"
collection: publications
category: manuscripts
permalink: /publication/F_test
excerpt: 'Рассматривается применение F-теста (критерия Фишера) в задаче регрессии'
date: 2026-04-10
venue: ''
slidesurl: 'https://academicpages.github.io/files/slides1.pdf'
paperurl: 'https://academicpages.github.io/files/paper1.pdf'
bibtexurl: ''
citation: ''
---
Рассматривается применение \\(F\\)-теста (критерия Фишера) в задаче регрессии. \\(F\\)-тест также широко используется для проверки равенства дисперсий двух или более выборок.

Предположим, что мы аппроксимируем данные полиномиальной функцией. Степень полинома, а следовательно, и количество параметров модели, можно увеличивать сколь угодно. Однако, следуя принципу бритвы Оккама, следует избегать чрезмерного усложнения модели. В противном случае возникает проблема переобучения: модель хорошо описывает обучающую выборку, но оказывается непригодной для предсказания новых данных. В связи с этим возникает задача проверки целесообразности усложнения модели.

Рассмотрим две модели, где модель 1 является вложенной в модель 2. Это означает, что модель 1 содержит \\(p_1\\) параметров, а модель 2 содержит \\(p_2\\) параметров, причём \\(p_1 < p_2\\). На тренировочных данных модель с большим числом параметров всегда может описать их не хуже, чем модель с меньшим числом параметров. Следовательно, модель 2, как правило, обеспечивает лучшее (то есть с меньшей ошибкой) соответствие данным. Чтобы определить, является ли это улучшение статистически значимым, используется \\(F\\)-тест.

Пусть имеется \\(n\\) наблюдений. Тогда \\(F\\)-статистика вычисляется по формуле:

\\[
F = \frac{(\chi^2_1 - \chi^2_2) / (p_2 - p_1)}{\chi^2_2 / (n - p_2)} = \frac{\chi^2_1 - \chi^2_2}{\chi^2_2} \cdot \frac{n - p_2}{p_2 - p_1},
\\]

где \\(\chi^2\\) — взвешенная сумма квадратов остатков, в которой веса равны обратным дисперсиям ошибок измерений.

При нулевой гипотезе \\(H_0\\): "модель 2 не обеспечивает лучшего соответствия данным по сравнению с моделью 1" - статистика \\(F\\) имеет \\(F\\)-распределение с \\((p_2 - p_1,\, n - p_2)\\) степенями свободы. Нулевая гипотеза отвергается, если наблюдаемое значение \\(F\\) превышает критическое значение при заданном уровне значимости.

```python
import numpy as np
from scipy import stats

# Generate synthetic dataset
np.random.seed(42)
x = np.linspace(-3, 3, 100)
y = 1 + 2*x + 0.5*x**2 + np.random.normal(0, 1, size=len(x))

alpha = 0.05          # significance level for F-test
max_degree = 5        # maximum polynomial degree to test
n_points = len(x)     # number of data points

chi_ls = []         # RSS or chi^2  

# Fit polynomial
def fit_rss(x, y, degree):
    # Build design matrix: [1, x, x^2, ..., x^degree]
    X = np.vander(x, degree + 1, increasing=True)
    # Solving
    beta = np.linalg.lstsq(X, y, rcond=None)[0]
    # Predicted values
    y_pred = X @ beta
    # Residual Sum of Squares
    rss = np.sum((y - y_pred)**2)
    return rss

# Compute RSS for models of increasing complexity
for deg in range(1, max_degree + 1):
    rss = fit_rss(x, y, deg)
    chi_ls.append(rss)

chi_ls = np.array(chi_ls)

# Number of parameters in each model
params = np.array([d + 1 for d in range(1, max_degree + 1)])

# Degrees of freedom
df = n_points - params

# Critical F value for each pair of nested models
f_critical = stats.f.ppf(1 - alpha, 1, df[1:])

# Compute F-statistics for consecutive model comparisons:
# (degree 1 vs 2), (2 vs 3), ...
f_stats = ((chi_ls[:-1] - chi_ls[1:]) / chi_ls[1:]) * (df[1:] / 1)

#Print results
print("chi^2 (RSS):", chi_ls)
for i, (f, f_cr) in enumerate(zip(f_stats, f_critical)):
    print(f"\nComparing degree {i+1} -> {i+2}")
    print(f"F = {f:.3f}, F_crit = {f_cr:.3f}")
    
    if f > f_cr:
        # Improvement is statistically significant → keep increasing complexity
        print("Model improvement is significant")
        best_degree = i + 2
    else:
        # No significant improvement → stop here
        print("No significant improvement, stop increasing degree")
        break

"""
# Output
chi^2 (RSS): [284.03466406  81.17986875  77.23541993  76.3083788   76.29129554]

Comparing degree 1 -> 2
F = 242.387, F_crit = 3.939
Model improvement is significant

Comparing degree 2 -> 3
F = 4.903, F_crit = 3.940
Model improvement is significant

Comparing degree 3 -> 4
F = 1.154, F_crit = 3.941
No significant improvement, stop increasing degree
"""
```

В моей работе по изучению спектральных и временных характеристик гамма-всплексков \\(F\\)-тест применялся для аппроксимации фона в кривых блеска, зарегистрированных прибором Fermi/GBM. Телескоп находится на околоземной орбите, что приводит к нестабильности фона во времени. В большинстве случаев достаточно аппроксимации нулевой или первой степени, однако иногда требуются модели более высокого порядка. Ниже приводятся несколько кривых блеска, которые имеют различные поведения фона, \\(F\\)-тест хорошо справляется с задачей определения оптимальной степени полинома без переусложнения модели:

<figure>
  <img src="/images/publications/SP1.png" alt="1 степень" style="height: auto;">
  <figcaption>Рис. 1: Аппроксимация полиномом 1-й степени</figcaption>
</figure>

<figure>
  <img src="/images/publications/SP2.png" alt="2 степень" style="height: auto;">
  <figcaption>Рис. 2: Аппроксимация полиномом 2-й степени</figcaption>
</figure>

<figure>
  <img src="/images/publications/SP3.png" alt="4 степень" style="height: auto;">
  <figcaption>Рис. 3: Аппроксимация полиномом 4-й степени</figcaption>
</figure>

