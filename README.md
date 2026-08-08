# FluidSound

Bubble-based water sound synthesis code based on the papers:
    
> [Improved Water Sound Synthesis using Coupled Bubbles](https://graphics.stanford.edu/papers/coupledbubbles/). Kangrui Xue, Ryan M. Aronson, Jui-Hsien Wang, Timothy R. Langlois, Doug L. James. *ACM Transactions on Graphics (SIGGRAPH North America 2023)*.  

> [Toward Animating Water with Complex Acoustic Bubbles](https://www.cs.cornell.edu/projects/Sound/bubbles/). Timothy R. Langlois, Changxi Zheng, Doug L. James. *ACM Transactions on Graphics (SIGGRAPH North America 2016)*. 

## Quick start: trackedBubInfo.txt → audio in one command

With `runFluidSound` built (see below), `generate_wav_from_trackedBubInfo.py` renders a
tracked-bubble file straight to a peak-normalized WAV **next to the input**, and can
optionally mux the audio onto a video. Two ready-to-run LBM scenes ship in `Scenes/`
(`Fruits` and `Exhale`: cleaned event graph, NN-predicted per-sample frequencies,
frequency x2 with matched radius /2). From this folder:

    # fruit splash -> Scenes/Fruits/trackedBubInfo_NN_fruit.{wav,mp4}
    python generate_wav_from_trackedBubInfo.py Scenes/Fruits/trackedBubInfo_NN_fruit.txt --video Scenes/Fruits/preview_3d_fruit_silent.mp4 --damping-coeff 0.8

    # underwater exhale -> Scenes/Exhale/trackedBubInfo_NN_exhale.{wav,mp4}
    python generate_wav_from_trackedBubInfo.py Scenes/Exhale/trackedBubInfo_NN_exhale.txt --video Scenes/Exhale/preview_3d_exhale_silent.mp4 --damping-coeff 0.8

Each scene also has a `trackedBubInfo_Minnaert_<scene>.txt` with the identical bubble
record and event graph but Minnaert frequencies — swap it in to A/B the frequency model.
`Scenes/GlassPour2016/` holds the WaveBlender reference pour from the papers.
Common knobs (defaults: scheme 0, 48 kHz, cutoff 10 ms, smoothstep, damping 1.0):

    python generate_wav_from_trackedBubInfo.py Scenes/Fruits/trackedBubInfo_NN_fruit.txt --scheme 1 --damping-coeff 0.7

## Quick start: separating N, S, M audio

`--forcing-types` renders per-event-type stems: only the listed bubble start-event
types (N = entrainment, S = split, M = merge) receive forcing impulses; all other
solver state (oscillator set, chaining, frequencies, damping, merge-forcing RNG
stream) is unchanged, so the three stems **sum sample-exactly to the full mix**:

    python generate_wav_from_trackedBubInfo.py Scenes/Fruits/trackedBubInfo_NN_fruit.txt --damping-coeff 0.8 --forcing-types N
    python generate_wav_from_trackedBubInfo.py Scenes/Fruits/trackedBubInfo_NN_fruit.txt --damping-coeff 0.8 --forcing-types S
    python generate_wav_from_trackedBubInfo.py Scenes/Fruits/trackedBubInfo_NN_fruit.txt --damping-coeff 0.8 --forcing-types M

(Each run overwrites `<stem>.wav` next to the input — pass `-o`-style separation by
copying the tracked file or moving the outputs between runs. Note each WAV is
peak-normalized individually; for loudness-true comparison rescale the stems by the
full mix's raw peak.)

`--event-log events.csv` additionally writes one row per forcing event with the
solver-measured per-impulse response:

    python generate_wav_from_trackedBubInfo.py Scenes/Fruits/trackedBubInfo_NN_fruit.txt --damping-coeff 0.8 --event-log Scenes/Fruits/events.csv

Columns: `t_event,osc_idx,bub_id,event_type,radius,cutoff_tau,weight,peak_accel` —
`weight`/`cutoff_tau` are the impulse parameters (0 weight = silenced by the split
mass guard, the merge two-parent condition, or `--forcing-types`), and `peak_accel`
is the maximum |v''| the integrator itself recorded while that impulse was the
oscillator's active one (exact from-rest birth response for a chain's first event;
includes inherited ringing for later chain links). The log is written after
synthesis completes.

The runFluidSound binary is auto-discovered (superproject `build/Release/`, an in-module
build, then PATH; `--exe` overrides), and `--video` finds ffmpeg via the imageio-ffmpeg
package or PATH (`--ffmpeg` overrides). Only numpy is required for the WAV itself. The
input must be reader-conformant (no `#` comment lines, resolvable Start/End references,
merge → 1 target, split → 2 children); frequency/radius scaling should be baked in
beforehand — see the full option list under *Command-line arguments* below.

## Build Instructions

**Dependencies:** C++11, Eigen 3.4

Building is handled by CMake. For example, to build from source on Mac & Linux:

    git clone https://github.com/kangruix/FluidSound
    cd FluidSound
    mkdir build && cd build && cmake -DCMAKE_BUILD_TYPE=Release ..
    make -j4

We provide an example scene in Scenes/GlassPour2016/. To run the code:

    ./runFluidSound ../Scenes/GlassPour2016/trackedBubInfo.txt 48000 1
    python ../scripts/write_wav.py output.txt 48000

With importance filtering (full graph + filtered subset for output):

    ./runFluidSound <full_bub_file> <filtered_bub_file> 48000 1 -o output_filtered.txt

Optional `-o`/`--output` sets the waveform output file (default: output.txt).

The full file provides the complete merge/split graph; the filtered file lists which bubbles contribute to the final sound. Bubbles not in the filtered set still participate in coupling but their contribution is excluded from the output. 

(where the 1 indicates the scheme: 0 - uncoupled, 1 - coupled). Afterwards, the simulated audio will be written to 'output.wav'. More scenes are available [here](https://graphics.stanford.edu/papers/waveblender/dataset/index.html).

## Command-line arguments

```
runFluidSound <bubFile> [filteredFile] <samplerate> <scheme> [options]
```

Positional arguments:

| Argument | Default | Description |
|---|---|---|
| `bubFile` | `../Scenes/GlassPour2016/trackedBubInfo.txt` | Full bubble graph (WaveBlender `trackedBubInfo` format; required for merge/split references). |
| `filteredFile` | *(none)* | Optional importance-filtered subset: only these bubbles contribute to the output, but all bubbles in `bubFile` still participate in coupling. |
| `samplerate` | `48000` | Audio sample rate in Hz. |
| `scheme` | `1` | `0` = uncoupled oscillators, `1` = coupled. |

Options:

| Option | Default | Description |
|---|---|---|
| `-o`, `--output FILE` | `output.txt` | Waveform output file (one sample per line; convert with `scripts/write_wav.py`). |
| `--max-time SEC` | off | If > 0, stop the integration at this time and write a partial waveform. |
| `--time-jitter W` | `0` | Shift each oscillator by U(−W, W) seconds to desynchronize grid clicks. |
| `--time-jitter-seed N` | `random_device` | RNG seed for `--time-jitter`. |
| `--transient-periods N` | `0` | Attenuate start forcing for oscillators lasting fewer than N periods. |
| `--transient-gain G` | `1` | Forcing multiplier applied to those transient oscillators. |
| `--forcing-cutoff SEC` | `0.0006` | Maximum duration of a start-forcing impulse. |
| `--forcing-envelope hard\|smoothstep` | `smoothstep` | Start impulse envelope shape. |
| `--dense-events` | **on** | Insert every per-sample-line solveTime into the integrator's event-time set, so K = ω₀² follows the trackedBubInfo frequency column at every row instead of a linear ramp from an oscillator's first to last sample. |
| `--no-dense-events` | — | Restore the legacy linear-ramp behavior. |
| `--damping-coeff C` | `1.0` | Multiplier on the per-sample β from `Oscillator::calcBeta` (unmodified Czerski/Deane radiative+viscous+thermal damping at 1.0). C < 1 lengthens the ringdown, C > 1 shortens it. |
| `--forcing-types STR` | `NSM` | Subset of `NSM`: only these bubble start-event types receive forcing impulses (N = entrain, S = split, M = merge). Impulses are zeroed after the forcing RNG draw, so stems rendered with `N` / `S` / `M` separately sum sample-exactly to the full `NSM` render. |
| `--event-log PATH` | off | Write a per-forcing-event CSV (`t_event,osc_idx,bub_id,event_type,radius,cutoff_tau,weight`) in final post-jitter/clip form; gated or guard-zeroed events appear with weight 0. |
| `--listener-position X Y Z` | off | Enable per-oscillator 1/distance attenuation, applied after integration (coupled dynamics untouched). Coordinates live in the bubble-position space of the input file. |
| `--listener-epsilon E` | `1e-6` | Lower clamp on the listener distance to avoid division blow-ups. |

All options also accept `snake_case` spellings (e.g. `--dense_events`).

### Miscellaneous

Documentation can be built using Doxygen: `doxygen Doxyfile`

### Citation

```bibtex
@article{Xue:2023:CoupledBubbles,
  author = {Xue, Kangrui and Aronson, Ryan M. and Wang, Jui-Hsien and Langlois, Timothy R. and James, Doug L.},
  title = {Improved Water Sound Synthesis using Coupled Bubbles},
  journal = {ACM Transactions on Graphics (Proceedings of SIGGRAPH 2023)},
  year = {2023}, volume = {42}, number = {4}, month = jul
}
@article{Langlois:2016:Bubbles,
  author = {Langlois, Timothy R. and Zheng, Changxi and James, Doug L.},
  title = {Toward Animating Water with Complex Acoustic Bubbles},
  journal = {ACM Transactions on Graphics (Proceedings of SIGGRAPH 2016)},
  year = {2016}, volume = {35}, number = {4}, month  = jul
}
```
