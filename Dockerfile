# Use the official Ubuntu image as the base
FROM ubuntu:latest

# Install any necessary dependencies
RUN apt-get update && apt-get install -y \
    git

# Clone the repository
RUN git clone --recurse-submodules https://github.com/cc2qe/simu.git /simu

# Set the working directory
WORKDIR /simu

# Install any additional dependencies specific to the software
RUN apt-get install -y \
    cmake \
    build-essential \
    libboost-all-dev

# (Optional) Set environment variables if needed
# ENV VAR_NAME=value

# (Optional) Run any setup scripts provided by the software
RUN bash ./setup.sh
RUN export PATH=$PATH:/simu/build/bin/simu

# (Optional) Expose any necessary ports
# EXPOSE <port>

# (Optional) Define the entry point or default command to run the software
# CMD ["<command>", "<arguments>"]

