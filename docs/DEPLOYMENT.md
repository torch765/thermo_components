# DigitalOcean App Platform Deployment

## Deployment Target

- Platform: DigitalOcean App Platform
- Source: GitHub repository `torch765/thermo_components`
- Branch: `main`
- Region: London
- Component type: Web Service
- Runtime: Python 3.13
- Initial size: 1 shared vCPU / 512 MiB
- Health check: `GET /health`

The repository-owned App Platform specification is
[`/.do/app.yaml`](../.do/app.yaml). The Python version is pinned in
[`/runtime.txt`](../runtime.txt).

## Initial Control Panel Setup

1. Open <https://cloud.digitalocean.com/apps>.
2. Select **Create App**.
3. Choose **GitHub** as the source provider.
4. Authorize DigitalOcean to read `torch765/thermo_components` if prompted.
5. Select:
   - Repository: `torch765/thermo_components`
   - Branch: `main`
   - Source directory: `/`
   - Autodeploy: enabled
6. Confirm that DigitalOcean detects one Python **Web Service**.
7. Review the service settings:
   - Name: `web`
   - Run command:

     ```text
     python -m uvicorn thermo_components.adapters.web.app:app --host 0.0.0.0 --port $PORT
     ```

   - HTTP port: `8080`
   - Health-check path: `/health`
8. Choose the London region.
9. Start with the 1 shared vCPU / 512 MiB plan.
10. Review the monthly price before selecting **Create Resources** or
    **Deploy**.

No database or environment variables are required for the MVP. App Platform
provides `PORT` automatically, and `lhv_data.db` is bundled read-only with the
repository.

## Resource Note

The initial service size is the current USD 5/month shared plan. The
thermodynamic dependency stack is heavier than a typical FastAPI application.
If runtime logs show an out-of-memory restart, change the service to the
1 shared vCPU / 1 GiB fixed plan before investigating application errors.

## Post-Deployment Smoke Test

After the deployment becomes active:

1. Open `/health` and confirm the JSON status is `ok`.
2. Open `/` and calculate 100% methane.
3. Open `/flow` and test `1 kg/h` to `t/d`.
4. Download the Excel report and open the workbook.
5. Review runtime logs for restarts, memory failures, or missing-resource
   errors.
