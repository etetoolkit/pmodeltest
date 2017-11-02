(for x in 2.7 3.4 3.5 3.6; do conda build --py $x ./meta.yaml; done)|grep "anaconda upload " > upload.sh
