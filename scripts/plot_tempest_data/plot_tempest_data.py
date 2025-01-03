import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import linregress, spearmanr
from sklearn.linear_model import LinearRegression


def create_model(input_file):

    df = pd.read_csv(input_file, sep='\t')

    X = df['date'].values.reshape(-1, 1)
    Y = df['distance'].values

    model = LinearRegression()
    model.fit(X, Y)

    Y_pred = model.predict(X)

    # Вычисление статистики линейной регрессии
    slope, intercept, r_value, p_value_linreg, std_err = linregress(df['date'], df['distance'])
    # Коэффициент корреляции Спирмена
    spearman_corr, p_value_spearman = spearmanr(df['date'], df['distance'])

    with open(rf'{work_dir}\{region_name}_p_value.csv', 'w') as pval:
        print('metric,values', file=pval)
        print(
            f"linear_regression,slope={slope};intercept={intercept};r-squared={r_value ** 2};p-value={p_value_linreg}",
            file=pval)
        print(f"spearman_correlation,coefficient={spearman_corr};p-value={p_value_spearman}",
              file=pval)

    create_plot(X, Y, Y_pred, input_file, p_value_linreg, p_value_spearman)


def create_plot(x, y, y_pred, input_file, p_val_linreg, p_val_spearman):

    plt.scatter(x, y, color='black', label='Data')
    plt.plot(x, y_pred, color='red', label='Linear Regression')
    plt.xlabel('Date')
    plt.ylabel('Root-to-tip divergence')
    plt.suptitle(f'   {region_name}', weight='bold')
    plt.title(f"\nLinear Regression p-value: {p_val_linreg:.3e}\nSpearman p-value: {p_val_spearman:.3e}",
              fontsize=7.5)
    plt.legend()

    plt.ylim(min(y) - 0.01, max(y) + 0.01)

    plt.savefig(rf'{work_dir}\{region_name}_root_to_tip_scatter.svg',
                format='svg')


work_dir = (r'C:\Users\gdvov\OneDrive\Documents\GitHub\astrovirus-main\trees\iqtree_output\tempest_analysis'
            r'\serotype_4\serotype_4_6048-7191')

region_name = work_dir.split('\\')[-1].split('.fas')[0]

create_model(rf'{work_dir}\serotype_4_6048-7191.fas.contree.tsv')
