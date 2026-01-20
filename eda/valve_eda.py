import marimo

__generated_with = "0.17.7"
app = marimo.App()


@app.cell
def _():
    import marimo as mo
    import pandas as pd
    return mo, pd


@app.cell
def _(pd):
    df = pd.read_csv('data/processed/reboot_times.csv')
    df
    return (df,)


@app.cell
def _(df, pd):
    # Plot the frequency of timestamp entries per day
    import matplotlib.pyplot as plt

    df['date'] = pd.to_datetime(df['timestamp']).dt.date
    frequency = df['date'].value_counts().sort_index()

    plt.figure(figsize=(10, 6))
    frequency.plot(kind='bar', color='skyblue')
    plt.title('Frequency of Timestamp Entries Per Day')
    plt.xlabel('Date')
    plt.ylabel('Number of Entries')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.show()
    return (frequency,)


@app.cell
def _(frequency):
    frequency
    return


@app.cell
def _(df):
    df
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## pump events

    1. Aug 18 141 reboots starting at 10:00 UTC, ending at 16:18 UTC
    2. Oct 1 41 reboots starting at 17:03, ending at 22:42
    """)
    return


if __name__ == "__main__":
    app.run()
