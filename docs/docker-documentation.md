# Using IBDMap with Docker
There is a Dockerfile in the root available that will build & prepare IBDMap to
run inside of a docker container based on Ubuntu 20.04. 

When built, the image can be used to run ibdmap by using 
[docker bind mounts](https://docs.docker.com/storage/bind-mounts/) to present
your files to the container. For example:

```shell
docker run -v $(pwd)/test:/app/test rjbohlender/ibdmap -i test/infile.txt.gz
```
will run ibdmap with the arguments "-i test/infile.txt.gz", and the contents of the 
local test directory mounted into the container at /app/test.

If you want to connect to the container interactively to debug:
```shell
docker run -it --entrypoint /bin/bash rjbohlender/ibdmap
```