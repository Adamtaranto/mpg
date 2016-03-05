from . import *

from docopt import docopt
from sys import stdout, stderr


def mpg_main():
    cli = '''
    USAGE:
        mpg [options] <reference> ...

    OPTIONS:
        -l LENGTH       Length of sequence to simulate [default: 0].
        -k ORDER        Markovian Order [default: 1].
        -s SEED         Random seed for generation. (uses /dev/urandom if unset)
        -d DUMPFILE     Dump data to DUMPFILE.
        -r              Reference file is a Yaml dump file.
    '''
    opts = docopt(cli)
    k = int(opts['-k'])
    l = int(opts['-l'])

    t = TransitionCounter(k)
    if opts['-r']:
        filename = opts['<reference>'][0]
        print('Loading reference dump from "{}"'.format(filename),
              file=stderr)
        t.load(filename)
    else:
        print('Inferring transition matrix from reference sequences',
              file=stderr)
        for fn in opts['<reference>']:
            t.consume_file(fn)
            print('Consumed', fn, file=stderr)

    if opts['-d']:
        filename = opts['-d']
        print('Saving reference dump to "{}"'.format(filename), file=stderr)
        t.save(filename)

    m = MarkovGenerator(t, seed=opts['-s'])
    print('Initialised Markov Generator', file=stderr)
    if l > 0:
        print('Generating sequence of {} bases'.format(l), file=stderr)
        print(seq2fa('genseq', m.generate_sequence(l)))
