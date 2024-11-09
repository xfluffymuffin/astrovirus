import pandas as pd
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression


def create_model(input_file):
    df = pd.read_csv(input_file, sep='\t')

    X = df['date'].values.reshape(-1, 1)
    Y = df['distance'].values

    model = LinearRegression()
    model.fit(X, Y)

    Y_pred = model.predict(X)

    create_plot(X, Y, Y_pred, input_file)


def create_plot(x, y, y_pred, input_file):
    region_name = input_file.split('\\')[-1].split('.fas')[0]

    plt.scatter(x, y, color='black', label='Data')
    plt.plot(x, y_pred, color='red', label='Linear Regression')
    plt.xlabel('Date')
    plt.ylabel('Root-to-tip divergence')
    plt.title(region_name)
    plt.legend()

    plt.ylim(min(y) - 0.01, max(y) + 0.01)

    plt.savefig(rf'saving_path\{region_name}_root_to_tip_scatter.svg',
                format='svg')


create_model(
    r"path_to_tsv_file_with_tempest_data\data.tsv")
