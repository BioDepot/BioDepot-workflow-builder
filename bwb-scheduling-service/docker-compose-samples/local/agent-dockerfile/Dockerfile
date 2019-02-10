FROM biodepot/bwb-scheduling-service
RUN apt-get update && apt-get install -y  \
    apt-transport-https ca-certificates curl software-properties-common gnupg2 \
    && curl -fsSL https://download.docker.com/linux/debian/gpg | apt-key add - \
    && apt-get remove -y curl gnupg2 \
    && apt-get autoclean -y \
    && apt-get autoremove -y \
    && rm -rf /var/lib/apt/lists/*

RUN add-apt-repository -y \
   "deb [arch=amd64] https://download.docker.com/linux/debian \
   stretch stable" \
   && apt-get update && apt-get install -y docker-ce docker-ce-cli \
   && apt-get remove -y apt-transport-https software-properties-common \
   && apt-get autoclean -y \
   && apt-get autoremove -y \
   && rm -rf /var/lib/apt/lists/*

