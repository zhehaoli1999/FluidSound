"""One-shot: trackedBubInfo.txt -> .wav, written next to the input file.

Runs runFluidSound on the given tracked file and converts the waveform to a
peak-normalized 16-bit WAV in the same directory:

    <input dir>/<input stem>_waveform.txt   (raw solver output, one sample/line)
    <input dir>/<input stem>.wav

The input must be reader-conformant (no '#' comment lines; every Start/End
reference resolvable; End M = 1 target, End S = 2 children) -- e.g. the output
of python/scene/_common/collapse_ghost_events_trackedbubinfo.py or of the
filter_scale_run pipeline. Frequencies and radii are used as-is: apply any
scaling beforehand (scale_tracked_bubinfo_frequency.py doubles frequency and
halves radius together).

Defaults follow the project's render settings rather than the bare solver
defaults: scheme 0 (uncoupled), 48 kHz, forcing cutoff 10 ms, smoothstep
envelope, dense events (solver default), damping coeff 1.0.

With --video, the rendered audio is additionally muxed onto the given mp4
(video stream copied, audio encoded as AAC, cut to the shorter of the two):

    <input dir>/<input stem>.mp4

Usage:
    python generate_wav_from_trackedBubInfo.py path/to/trackedBubInfo.txt
    python generate_wav_from_trackedBubInfo.py tracked.txt --scheme 1 --damping-coeff 0.7
    python generate_wav_from_trackedBubInfo.py tracked.txt --video preview_3d.mp4
    python generate_wav_from_trackedBubInfo.py tracked.txt --exe /path/to/runFluidSound

Only numpy is required (WAV writing uses the standard library); --video needs
an ffmpeg binary (the imageio-ffmpeg bundled one is found automatically if the
package is installed, else ffmpeg from PATH, else pass --ffmpeg).
"""

from __future__ import annotations

import argparse
import subprocess
import sys
import wave
from pathlib import Path

import numpy as np

_THIS_DIR = Path(__file__).resolve().parent

# Candidate runFluidSound binaries, most specific first: the BubbleGym
# superproject build, then an in-module build, then whatever is on PATH.
_EXE_CANDIDATES = [
    _THIS_DIR.parent.parent / "build" / "Release" / "runFluidSound.exe",
    _THIS_DIR.parent.parent / "build" / "runFluidSound",
    _THIS_DIR / "build" / "Release" / "runFluidSound.exe",
    _THIS_DIR / "build" / "runFluidSound",
    Path("runFluidSound.exe"),
    Path("runFluidSound"),
]


def find_exe(override: Path | None) -> Path:
    if override is not None:
        if not override.is_file():
            raise SystemExit(f"--exe not found: {override}")
        return override
    for cand in _EXE_CANDIDATES:
        if cand.is_file():
            return cand
    raise SystemExit(
        "runFluidSound executable not found. Build it (see README.md) or pass --exe.\n"
        "Tried:\n  " + "\n  ".join(str(c) for c in _EXE_CANDIDATES)
    )


def find_ffmpeg(override: Path | None) -> str:
    if override is not None:
        if not override.is_file():
            raise SystemExit(f"--ffmpeg not found: {override}")
        return str(override)
    try:
        import imageio_ffmpeg  # type: ignore

        return imageio_ffmpeg.get_ffmpeg_exe()
    except Exception:
        pass
    import shutil

    exe = shutil.which("ffmpeg")
    if exe:
        return exe
    raise SystemExit(
        "ffmpeg not found for --video muxing: pip install imageio-ffmpeg, "
        "put ffmpeg on PATH, or pass --ffmpeg."
    )


def mux(video: Path, wav_path: Path, out_mp4: Path, ffmpeg: str) -> None:
    cmd = [ffmpeg, "-y", "-loglevel", "error",
           "-i", str(video), "-i", str(wav_path),
           "-c:v", "copy", "-c:a", "aac", "-shortest", str(out_mp4)]
    print("+", " ".join(cmd))
    result = subprocess.run(cmd)
    if result.returncode != 0 or not out_mp4.is_file():
        raise SystemExit(f"ffmpeg mux failed with exit code {result.returncode}")
    print(f"muxed -> {out_mp4}")


def write_wav(waveform_txt: Path, wav_path: Path, samplerate: int) -> None:
    data = np.loadtxt(waveform_txt)
    if data.ndim != 1 or data.size == 0:
        raise SystemExit(f"unexpected waveform contents in {waveform_txt}")
    peak = float(np.max(np.abs(data)))
    if peak > 0:
        data = data / peak
    pcm = np.clip(data * 32767.0, -32768, 32767).astype("<i2")
    with wave.open(str(wav_path), "wb") as w:
        w.setnchannels(1)
        w.setsampwidth(2)
        w.setframerate(samplerate)
        w.writeframes(pcm.tobytes())
    print(f"{data.size / samplerate:.2f} s | peak {peak:.4g} (normalized) -> {wav_path}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    ap.add_argument("tracked", type=Path, help="trackedBubInfo-style input file")
    ap.add_argument("--exe", type=Path, default=None, help="runFluidSound executable override")
    ap.add_argument("--samplerate", type=int, default=48000)
    ap.add_argument("--scheme", type=int, default=0, choices=[0, 1],
                    help="0 = uncoupled (default here), 1 = coupled")
    ap.add_argument("--forcing-cutoff", type=float, default=0.01,
                    help="max start-impulse duration in seconds (default 0.01)")
    ap.add_argument("--forcing-envelope", choices=["smoothstep", "hard"], default="smoothstep")
    ap.add_argument("--damping-coeff", type=float, default=1.0)
    ap.add_argument("--no-dense-events", action="store_true",
                    help="legacy linear frequency ramp per oscillator (dense events are on by default)")
    ap.add_argument("--max-time", type=float, default=None,
                    help="stop the integration at this time (seconds)")
    ap.add_argument("--video", type=Path, default=None,
                    help="optional mp4: mux the rendered audio onto it, writing "
                         "<input stem>.mp4 next to the tracked input")
    ap.add_argument("--ffmpeg", type=Path, default=None,
                    help="ffmpeg executable override (used only with --video)")
    args = ap.parse_args()

    tracked = args.tracked.resolve()
    if not tracked.is_file():
        raise SystemExit(f"input not found: {tracked}")
    if args.video is not None and not args.video.is_file():
        raise SystemExit(f"--video not found: {args.video}")
    exe = find_exe(args.exe)
    ffmpeg = find_ffmpeg(args.ffmpeg) if args.video is not None else None

    waveform_txt = tracked.parent / f"{tracked.stem}_waveform.txt"
    wav_path = tracked.parent / f"{tracked.stem}.wav"

    cmd = [str(exe), str(tracked), str(args.samplerate), str(args.scheme),
           "-o", str(waveform_txt),
           "--forcing-cutoff", str(args.forcing_cutoff),
           "--forcing-envelope", args.forcing_envelope,
           "--damping-coeff", str(args.damping_coeff)]
    if args.no_dense_events:
        cmd.append("--no-dense-events")
    if args.max_time is not None:
        cmd.extend(["--max-time", str(args.max_time)])

    print("+", " ".join(cmd))
    result = subprocess.run(cmd)
    if result.returncode != 0:
        raise SystemExit(f"runFluidSound failed with exit code {result.returncode}")

    write_wav(waveform_txt, wav_path, args.samplerate)
    if args.video is not None:
        assert ffmpeg is not None
        mux(args.video.resolve(), wav_path, tracked.parent / f"{tracked.stem}.mp4", ffmpeg)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
