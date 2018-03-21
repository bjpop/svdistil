#!/bin/bash

set -e
errors=0

# Run unit tests
python svdistil/svdistil_test.py || {
    echo "'python python/svdistil/svdistil_test.py' failed"
    let errors+=1
}

# Check program style
pylint -E svdistil/*.py || {
    echo 'pylint -E svdistil/*.py failed'
    let errors+=1
}

[ "$errors" -gt 0 ] && {
    echo "There were $errors errors found"
    exit 1
}

echo "Ok : Python specific tests"
