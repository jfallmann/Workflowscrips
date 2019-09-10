# AddStructure.py ---
#
# Filename: AddStructure.py
# Description:
# Author: Joerg Fallmann
# Maintainer:
# Created: Tue Sep 10 18:00:42 2019 (+0200)
# Version:
# Package-Requires: ()
# Last-Updated: Tue Sep 10 18:02:18 2019 (+0200)
#           By: Joerg Fallmann
#     Update #: 1
# URL:
# Doc URL:
# Keywords:
# Compatibility:
#
#

# Commentary:
#
#
#
#

# Change Log:
#
#
#
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or (at
# your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with GNU Emacs.  If not, see <http://www.gnu.org/licenses/>.
#
#

# Code:

import sys,os
import RNA

for line in sys.stdin:
    sequence = '\t'.split(line)[-1]
    # create new fold_compound object
    fc = RNA.fold_compound(sequence)
    # compute minimum free energy (mfe) and corresponding structure
    (ss, mfe) = fc.mfe()

    sys.stdout.write('\t'.join([line.strip(),ss])+'\n')
#
# AddStructure.py ends here
