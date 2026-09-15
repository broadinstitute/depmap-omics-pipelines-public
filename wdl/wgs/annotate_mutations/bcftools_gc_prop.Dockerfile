FROM rust:1.81-bullseye AS gc_prop_builder

RUN apt-get -y update \
    && apt-get -y dist-upgrade \
    && apt-get install -y --no-install-recommends --no-install-suggests libclang-dev

WORKDIR /usr/src

RUN cargo new --bin gc_prop_in_window
WORKDIR /usr/src/gc_prop_in_window

COPY ./gc_prop_in_window/Cargo.toml ./gc_prop_in_window/Cargo.lock ./
RUN cargo build --release
RUN rm -rf src

COPY ./gc_prop_in_window/src ./src
RUN rm ./target/release/deps/gc_prop_in_window*
RUN cargo build --release

FROM us-central1-docker.pkg.dev/depmap-omics/terra-images/bcftools:production

WORKDIR /app
COPY --from=gc_prop_builder /usr/src/gc_prop_in_window/target/release/gc_prop_in_window .
