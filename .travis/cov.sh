cd ..

mkdir cov

cp src/test/CMakeFiles/*.dir/*.gc?? cov

cd cov

bash <(curl -s https://codecov.io/bash)