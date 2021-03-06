import subprocess
import logging
logging.getLogger().setLevel(logging.INFO)

import xarray as xr
import matplotlib.pyplot as plt
import cmocean
import ffmpeg

SMALL_SIZE = 12
MEDIUM_SIZE = 14
BIGGER_SIZE = 16

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

def make_movie(filename_pattern, movie_filename, fps=15):
    (
        ffmpeg
        .input(filename_pattern, framerate=fps)
        .output(movie_filename, crf=15, pix_fmt='yuv420p')
        .overwrite_output()
        .run()
    )

def plot_forcing(ds, filename="forcing_time_series.png"):
    fig, axes = plt.subplots(nrows=3, figsize=(16, 9), dpi=300)
    lat, lon = ds.attrs["lat"], ds.attrs["lon"]
    fig.suptitle(f"Forcing at {lat}°N, {lon}°E")

    time = ds.time.values / 86400
    τx = ds.τx.values
    τy = ds.τy.values
    QT = ds.QT.values
    QS = ds.QS.values

    axes[0].plot(time, τx, label="τx")
    axes[0].plot(time, τy, label="τy")
    axes[0].legend(loc="best", frameon=False)
    axes[0].set_xlabel("")
    axes[0].set_ylabel("Wind stress ()N/m²)")

    axes[1].plot(time, QT)
    axes[1].set_xlabel("")
    axes[1].set_ylabel("QT ()W/m²)")

    axes[2].plot(time, QS)
    axes[2].set_xlabel("time (days)")
    axes[2].set_ylabel("QS (kg/m²/s)")

    for ax in axes:
        ax.label_outer()

    logging.info(f"Saving: {filename}...")
    plt.savefig(filename)

def plot_large_scale(ds):
    for n in range(0, ds.time.size, 10):
        fig, axes = plt.subplots(ncols=3, figsize=(16, 9), dpi=300)

        t = ds.time.values[n] / 86400
        lat, lon = ds.attrs["lat"], ds.attrs["lon"]
        fig.suptitle(f"Large scale at {lat}°N, {lon}°E, t = {t:.3f} days", fontsize=16)

        Ugeo = ds.Ugeo.isel(time=n).squeeze()
        Vgeo = ds.Vgeo.isel(time=n).squeeze()
        Ugeo.plot(ax=axes[0], y="zC", label="Ugeo", color="tab:green")
        Vgeo.plot(ax=axes[0], y="zC", label="Vgeo", color="tab:red")

        u = ds.u.isel(time=n).squeeze()
        v = ds.v.isel(time=n).squeeze()
        u.plot(ax=axes[0], y="zC", label="U", color="tab:green", linestyle="--")
        v.plot(ax=axes[0], y="zC", label="V", color="tab:red", linestyle="--")

        axes[0].legend(loc="best", frameon=False)
        axes[0].set_xlabel("velocity [m/s]")

        T = ds.T.isel(time=n).squeeze()
        T.plot(ax=axes[1], y="zC")

        S = ds.S.isel(time=n).squeeze()
        S.plot(ax=axes[2], y="zC")

        for ax in axes:
            ax.set_title("")
            ax.label_outer()

        png_filename = f"large_scale_{n//10:05d}.png"
        logging.info(f"Saving: {png_filename}...")
        plt.savefig(png_filename)

        plt.close("all")

def plot_statistics(ds):
    for n in range(0, ds.time.size):
        fig, axes = plt.subplots(nrows=2, ncols=4, sharey=True, figsize=(16, 9), dpi=300)

        t = ds.time.values[n] / 86400
        lat, lon = ds.attrs["lat"], ds.attrs["lon"]
        fig.suptitle(f"Statistics at {lat}°N, {lon}°E, t = {t:.3f} days", fontsize=16)

        u = ds.u.isel(time=n).squeeze()
        v = ds.v.isel(time=n).squeeze()
        u.plot(ax=axes[0, 0], y="zC", label="u", color="tab:green")
        v.plot(ax=axes[0, 0], y="zC", label="v", color="tab:red")
        axes[0, 0].legend(loc="best", frameon=False)
        axes[0, 0].set_xlabel("velocity [m/s]")

        T = ds.T.isel(time=n).squeeze()
        T.plot(ax=axes[0, 1], y="zC")

        S = ds.S.isel(time=n).squeeze()
        S.plot(ax=axes[0, 2], y="zC")

        uu = ds.uu.isel(time=n).squeeze()
        vv = ds.vv.isel(time=n).squeeze()
        ww = ds.ww.isel(time=n).squeeze()
        uu.plot(ax=axes[1, 0], y="zC", label="$\overline{u^\prime u^\prime}$")
        vv.plot(ax=axes[1, 0], y="zC", label="$\overline{v^\prime v^\prime}$")
        ww.plot(ax=axes[1, 0], y="zC", label="$\overline{w^\prime w^\prime}$")
        axes[1, 0].legend(loc="best", frameon=False)
        axes[1, 0].set_xlabel("turbulent kinetic energy [m²/s²]")

        uv = ds.uv.isel(time=n).squeeze()
        uw = ds.uw.isel(time=n).squeeze()
        vw = ds.vw.isel(time=n).squeeze()
        uv.plot(ax=axes[1, 1], y="zC", label="$\overline{u^\prime v^\prime}$")
        uw.plot(ax=axes[1, 1], y="zC", label="$\overline{u^\prime w^\prime}$")
        vw.plot(ax=axes[1, 1], y="zC", label="$\overline{v^\prime w^\prime}$")
        axes[1, 1].legend(loc="best", frameon=False)
        axes[1, 1].set_xlabel("velocity covariance [m²/s²]")

        wT = ds.wT.isel(time=n).squeeze()
        wT.plot(ax=axes[1, 2], y="zC")
        axes[1, 2].set_xlabel("$\overline{w^\prime T^\prime}$")

        wS = ds.wS.isel(time=n).squeeze()
        wS.plot(ax=axes[1, 3], y="zC")
        axes[1, 3].set_xlabel("$\overline{w^\prime S^\prime}$")

        for ax in axes.reshape(-1):
            ax.set_title("")

        png_filename = f"statistics_{n:05d}.png"
        logging.info(f"Saving: {png_filename}...")
        plt.savefig(png_filename)

        plt.close("all")

def plot_slices(ds):
    for n in range(0, ds.time.size):
        fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(16, 9), dpi=300)

        t = ds.time.values[n] / 86400
        lat, lon = ds.attrs["lat"], ds.attrs["lon"]
        fig.suptitle(f"Slices at {lat}°N, {lon}°E, t = {t:.3f} days", fontsize=16)

        u = ds.u.isel(time=n).squeeze()
        u.plot.pcolormesh(ax=axes[0, 0], vmin=-0.5, vmax=0.5, cmap=cmocean.cm.balance, extend="both")
        axes[0, 0].set_title("")

        w = ds.w.isel(time=n).squeeze()
        w.plot.pcolormesh(ax=axes[0, 1], vmin=-0.5, vmax=0.5, cmap=cmocean.cm.balance, extend="both")
        axes[0, 1].set_title("")

        T = ds.T.isel(time=n).squeeze()
        T.plot.pcolormesh(ax=axes[1, 0], cmap=cmocean.cm.thermal, extend="both")
        axes[1, 0].set_title("")

        S = ds.S.isel(time=n).squeeze()
        S.plot.pcolormesh(ax=axes[1, 1], cmap=cmocean.cm.haline, extend="both")
        axes[1, 1].set_title("")

        png_filename = f"slice_{n:05d}.png"
        logging.info(f"Saving: {png_filename}...")
        plt.savefig(png_filename)

        plt.close("all")

def plot_surfaces(ds):
    for n in range(0, ds.time.size):
        fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(16, 9), dpi=300)

        t = ds.time.values[n] / 86400
        lat, lon = ds.attrs["lat"], ds.attrs["lon"]
        fig.suptitle(f"Surface at {lat}°N, {lon}°E, t = {t:.3f} days", fontsize=16)

        u = ds.u.isel(time=n).squeeze()
        u.plot.pcolormesh(ax=axes[0, 0], vmin=-0.5, vmax=0.5, cmap=cmocean.cm.balance, extend="both")
        axes[0, 0].set_title("")

        w = ds.w.isel(time=n).squeeze()
        w.plot.pcolormesh(ax=axes[0, 1], vmin=-0.5, vmax=0.5, cmap=cmocean.cm.balance, extend="both")
        axes[0, 1].set_title("")

        T = ds.T.isel(time=n).squeeze()
        T.plot.pcolormesh(ax=axes[1, 0], cmap=cmocean.cm.thermal, extend="both")
        axes[1, 0].set_title("")

        S = ds.S.isel(time=n).squeeze()
        S.plot.pcolormesh(ax=axes[1, 1], cmap=cmocean.cm.haline, extend="both")
        axes[1, 1].set_title("")

        png_filename = f"surface_{n:05d}.png"
        logging.info(f"Saving: {png_filename}...")
        plt.savefig(png_filename)

        plt.close("all")

def make_lesbrary_plots(lat, lon, days):
    dsf = xr.open_dataset(f"lesbrary_lat{lat}_lon{lon}_days{days}_fields.nc")
    dss = xr.open_dataset(f"lesbrary_lat{lat}_lon{lon}_days{days}_surface.nc")
    dsx = xr.open_dataset(f"lesbrary_lat{lat}_lon{lon}_days{days}_slice.nc")
    dsp = xr.open_dataset(f"lesbrary_lat{lat}_lon{lon}_days{days}_profiles.nc")
    dsl = xr.open_dataset(f"lesbrary_lat{lat}_lon{lon}_days{days}_large_scale.nc")

    plot_forcing(dsl)

    plot_large_scale(dsl)
    make_movie("large_scale_%05d.png", "large_scale.mp4")

    plot_statistics(dsp)
    make_movie("statistics_%05d.png", "statistics.mp4")

    plot_slices(dsx)
    make_movie("slice_%05d.png", "slice.mp4")

    plot_surfaces(dss)
    make_movie("surface_%05d.png", "surface.mp4")

make_lesbrary_plots(-50, 275, 10)

