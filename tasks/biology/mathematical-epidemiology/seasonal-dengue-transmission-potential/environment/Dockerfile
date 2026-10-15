FROM python:3.12.8-slim-bookworm

ENV PYTHONDONTWRITEBYTECODE=1 \
    PIP_NO_CACHE_DIR=1 \
    OMP_NUM_THREADS=1 \
    OPENBLAS_NUM_THREADS=1

# All packages, including transitive dependencies, are pinned.
RUN pip install --no-cache-dir --no-deps \
        numpy==2.2.5 \
        scipy==1.15.3 \
        pandas==2.2.3 \
        python-dateutil==2.9.0.post0 \
        pytz==2025.2 \
        tzdata==2025.2 \
        six==1.17.0 \
        pytest==8.4.1 \
        iniconfig==2.1.0 \
        packaging==24.2 \
        pluggy==1.6.0 \
        pygments==2.19.1 \
    && pip check

WORKDIR /app
COPY data/ /app/data/
RUN mkdir -p /app/output
