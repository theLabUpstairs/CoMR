# Building CoMR Container Images

This guide is intended for maintainers who need to produce the Docker image
(`comr:latest` by default) and the Singularity/Apptainer `.sif` artifact that
end users will run. If you simply want to execute CoMR using prebuilt images,
see `README.md`.

## Prerequisites

- Linux host with Docker (or compatible OCI runtime) and outbound HTTP/FTP
  access.
- Ability to download the licensed TargetP archive separately (not bundled in
  the public image; users must provide it at run time).
- Optional: Singularity/Apptainer installation if you plan to build `.sif`
  locally instead of using Sylabs Cloud.

## 1. Fetch third-party assets

Run the helper script once to download and unpack MitoFates and Mitoprot II:

```bash
bash scripts/fetch_third_party.sh
```

This populates `third_party/MitoFates` and `third_party/mitoprotII`, which are
copied into `/opt/software` inside the image build.

## 2. Build the Docker image

From the repository root:

```bash
docker build -t comr:latest .
```

If outbound network access is restricted, download the tarballs referenced in
`Dockerfile` ahead of time and adjust the `COPY`/`ADD` paths accordingly.

## 3. Build a Singularity/Apptainer image

Most HPC clusters require a `.sif` and disallow local Docker builds. 

```bash
docker build -t comr:latest .
singularity build CoMR.sif docker-daemon://comr:latest
```

This requires Singularity/Apptainer to be installed on the same machine that
already has Docker (or another daemon exposing the `docker-daemon://` API).

Remember that TargetP remains user-supplied at runtime due to licensing, so do
not bake it into the public image. 
