# JCORRAN Library
Contact: Maxim Virta <maxim.virta@cern.ch>, DongJo Kim <dongjo.kim@cern.ch>

Jyvaskyla Correlation analysis library package for Linux and macOS systems.

Recent papers to produce the results:
- Phys.Lett.B 835 (2022) 137485: New constraints for QCD matter from improved Bayesian parameter estimation in heavy-ion collisions at LHC
- arXiv:2411.01932: Enhancing Bayesian parameter estimation by adapting to multiple energy scales in RHIC and LHC heavy-ion collisions

## Installation

1. Clone this repository
2. Add the "bin"-folder to your PATH environment variable (see below)
3. Run `createJCORRANlib` to build the library

### Adding bin-folder to PATH

For Linux and macOS, add the following line to your shell profile file (`.bashrc`, `.zshrc`, etc.):

```bash
export PATH="/absolute/path/to/JCORRAN/bin:$PATH"
```

Then reload your profile with: 

For bash:
```bash
source ~/.bashrc
```

For zsh:
```bash
source ~/.zshrc
```

# Commands
`jcorran-config-path` outputs the path to JCORRAN folder.

`jcorran-config-lib` outputs the lib flags required for linking.

`jcorran-config-inc` outputs the cxx flags required for compiling.