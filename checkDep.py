import warnings
import re


def ver2int(v):
    return int("".join(v.split('.')))

not_installed = []

for lname, lver in [('numpy','1.11.0'),('pandas','0.17.1'),('matplotlib','1.5.1')]:
    try:
        mod = __import__(lname)
        if ver2int(mod.__version__) < ver2int(lver):
            warnings.warn('Installed {0} version ({1}) is older than {2}! Some functions may not work as expected...'.format(lname,mod.__version__,lver))
    except ImportError:
        not_installed.append(lname)

try:
    import intervaltree as itree
    itree_version = re.findall(r'Version\s*([\d.]+)',itree.__doc__)[0]
    if ver2int(itree_version) < ver2int('2.0'):
        warnings.warn('Installed intervaltree version ({0}) is older than 2.0! Some functions may not work as expected...'.format(itree_version))
except ImportError:
    not_installed.append('intervaltree (https://pypi.python.org/pypi/intervaltree)')

if not_installed:
    print "Following libraries are not installed:"
    for ni in not_installed:
        print ni
    print "to use this repo first install all libraries listed above!"
else:
    print "All libraries installed!"
