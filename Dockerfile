# Lightweight Docker container for Boltzgen
# Provides a standardized Linux environment for testing and inference with torch.compile support

# CUDA / cuDNN base with no Python
FROM nvidia/cuda:12.8.0-cudnn-devel-ubuntu24.04

# System prerequisites + Python 3.12
# Note: Modal uses Python 3.10, but we use 3.12 for better compatibility
ENV DEBIAN_FRONTEND=noninteractive \
    PYTHON_VERSION=3.12.7 \
    PATH=/usr/local/bin:$PATH \
    TF_CPP_MIN_LOG_LEVEL=2 \
    TF_ENABLE_ONEDNN_OPTS=0 \
    TOKENIZERS_PARALLELISM=true

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        build-essential curl git ca-certificates wget \
        libssl-dev zlib1g-dev libbz2-dev libreadline-dev \
        libsqlite3-dev libncursesw5-dev xz-utils tk-dev \
        libxml2-dev libxmlsec1-dev libffi-dev liblzma-dev \
        ninja-build && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

RUN curl -fsSLO https://www.python.org/ftp/python/${PYTHON_VERSION}/Python-${PYTHON_VERSION}.tgz && \
    tar -xzf Python-${PYTHON_VERSION}.tgz && \
    cd Python-${PYTHON_VERSION} && \
    ./configure --enable-optimizations && \
    make -j"$(nproc)" && \
    make altinstall && \
    cd .. && rm -rf Python-${PYTHON_VERSION}* && \
    ln -s /usr/local/bin/python3.12 /usr/local/bin/python && \
    ln -s /usr/local/bin/pip3.12 /usr/local/bin/pip

# Location of project code (inside image) – NOT shared with host
WORKDIR /app

# Copy requirements first for layer caching (matching Modal image order)
COPY requirements.txt .
COPY requirements_modal.txt .

# Install packages
# Order is important
# 1. Upgrade pip/setuptools
# 2. Install this repo
# 3. Install requirements.txt
# 4. Force reinstall torch/torchvision with CUDA 12.8, so we have the most up to date version
# 5. Force reinstall numpy, numpy > 2.0 fails with scipy and some other packages, so we manually revert to 1.26.4
RUN pip install --upgrade pip setuptools && \
    pip install --no-cache-dir -e /app && \
    pip install -r requirements.txt && \
    pip install --force-reinstall torch torchvision --index-url https://download.pytorch.org/whl/cu128 -U && \
    pip install --force-reinstall numpy==1.26.4

# Copy the rest of the source
COPY . .

# ──────────────────────────────────────────────────────────────────────────────
# Single persistent host volume (/workspace) for *all* artefacts & caches
# Bind-mount it when you run the container:  -v ${PWD}:/workspace
# ──────────────────────────────────────────────────────────────────────────────
ENV PROJECT_ROOT=/workspace \
    PYTHONPATH=/app \
    CHAI_DOWNLOADS_DIR=/workspace/models/chai1 \
    HF_HUB_ENABLE_HF_TRANSFER=1 \
    DISABLE_PANDERA_IMPORT_WARNING=True \
    HF_HOME=/workspace/.cache/huggingface \
    TORCH_HOME=/workspace/.cache/torch \
    XDG_CACHE_HOME=/workspace/.cache \
    WANDB_DIR=/workspace/logs \
    TQDM_CACHE=/workspace/.cache/tqdm

RUN mkdir -p \
      /workspace/.cache/huggingface \
      /workspace/.cache/torch \
      /workspace/.cache/tqdm \
      /workspace/logs \
      /workspace/data \
      /workspace/results

ARG DOWNLOAD_WEIGHTS=false
RUN mkdir -p "${HF_HOME}" && \
    if [ "${DOWNLOAD_WEIGHTS}" = "true" ]; then \
        boltzgen download all --cache "${HF_HOME}" --force_download; \
    fi

ARG USERNAME=boltzgen
ARG USER_UID=1000
ARG USER_GID=1000

RUN groupadd --gid ${USER_GID} ${USERNAME} && \
    useradd --uid ${USER_UID} --gid ${USER_GID} --create-home --shell /bin/bash ${USERNAME}

RUN mkdir -p "${HF_HOME}" && chown -R ${USER_UID}:${USER_GID} "${HF_HOME}"

USER ${USERNAME}
WORKDIR /workspace

# Declare the volume so other developers know it's intended to persist
VOLUME ["/workspace"]

ENTRYPOINT ["boltzgen"]
CMD ["--help"]
