#!/usr/bin/env python3
"""Compare (or set) the simulator's compiled flag defaults against the
production argv in cpp/production/prodargs.sh.

Production flags must equal the compiled defaults, except the per-run values
the caller sets (seed, output file, step counts).  Shell variables in prodargs
resolve to the production defaults below.

usage: python sync_defaults.py [--apply] [--prodargs PATH] [--main PATH]
Exit status 1 when a production flag differs from its compiled default.
"""
import argparse
import re
import sys
from pathlib import Path

here = Path(__file__).resolve().parent
PER_RUN = {'seed', 'filename', 'nsteps', 'save_every'}
# Values of the shell variables in prodargs.sh for the production campaign.
SHELL_DEFAULTS = {'HPH': '27', 'KILL': '0.050', 'RATE': '1', 'SLAB': 'false', 'CUTOFF': '0'}
OPTION = re.compile(
    r'\("(?P<name>[A-Za-z_0-9]+)",\s*"(?:[^"\\]|\\.)*",\s*cxxopts::value<[^>]+>\([^)]*\)'
    r'(?P<chain>(?:->[a-z_]+\("[^"]*"\))*)')


def production_flags(path):
    flags = {}
    for name, value in re.findall(r'--([A-Za-z_0-9]+)=(\S+)', path.read_text()):
        variable = re.fullmatch(r'"?\$\{(\w+)\}"?', value)
        flags[name] = SHELL_DEFAULTS.get(variable.group(1)) if variable else value
    return flags


def same(a, b):
    try:
        return abs(float(a) - float(b)) <= 1e-12 * max(1.0, abs(float(a)))
    except ValueError:
        return a.strip('"').lower() == b.strip('"').lower()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--apply', action='store_true', help='rewrite main.cpp defaults')
    parser.add_argument('--prodargs', type=Path, default=here/'prodargs.sh')
    parser.add_argument('--main', type=Path, default=here.parent/'src'/'main.cpp')
    args = parser.parse_args()
    flags = production_flags(args.prodargs)
    source = args.main.read_text()
    defaults = {}
    for match in OPTION.finditer(source):
        found = re.search(r'default_value\("([^"]*)"\)', match.group('chain'))
        defaults[match.group('name')] = found.group(1) if found else None
    unknown = sorted(set(flags) - set(defaults))
    differing = [(name, value) for name, value in flags.items()
                 if name not in PER_RUN and name in defaults and value is not None
                 and (defaults[name] is None or not same(value, defaults[name]))]
    for name, value in differing:
        print(f'{name}: production {value}, compiled default {defaults[name]}')
    if unknown:
        print('unknown to main.cpp:', ', '.join(unknown))
    if not args.apply:
        print(f'{len(differing)} production flags differ from the compiled defaults')
        return 1 if differing or unknown else 0

    def rewrite(match):
        name = match.group('name')
        if name in PER_RUN or name not in flags or flags[name] is None:
            return match.group(0)
        value = flags[name].strip('"')
        text = match.group(0)
        assert 'default_value("' in text, f'{name} has no default_value'
        return re.sub(r'default_value\("[^"]*"\)', f'default_value("{value}")', text, count=1)

    args.main.write_text(OPTION.sub(rewrite, source))
    print(f'set {len(differing)} compiled defaults to production values in {args.main}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
