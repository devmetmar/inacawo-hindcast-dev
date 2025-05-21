class getERA5:
    def __init__(self, tstart:str, tend:str):
        self.tstart = tstart
        self.tend = tend
        self.logger = self._setup_logger()
        self.logger.info(f"""
        ============================
        Initialized ERA5 downloader.
        TSTART: {self.tstart}
        TEND  : {self.tend}
        ============================
        """)
        self.config = self.load_yaml()

    def _setup_logger(self):
        import os
        import logging
        from datetime import datetime
        
        logdir  = datetime.strptime(self.tstart, "%Y%m%d").strftime("/home/cawofcst_dev/devmetmar/logs/download/era5/%Y/%m")
        logname = f"{logdir}/era5_{self.tstart}_{self.tend}_{datetime.utcnow().strftime('%Y%m%d_%H%M%S')}.log"
        if not os.path.exists(logdir):
            os.makedirs(logdir)
        # logging.basicConfig(
        #     filename=logname,
        #     level=logging.INFO,
        #     format="%(asctime)s - %(levelname)s - %(message)s"
        # )
        logger = logging.getLogger("ERA5Downloader")
        logger.setLevel(logging.INFO)

        if not logger.handlers:
            # File handler
            file_handler = logging.FileHandler(logname)
            file_handler.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))
            logger.addHandler(file_handler)

            # Console handler
            console_handler = logging.StreamHandler()
            console_handler.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))
            logger.addHandler(console_handler)

        return logger

    def load_yaml(self, config="/home/cawofcst_dev/devmetmar/cawo-anal-202412/scripts/config.yml"):
        import yaml
        from pprint import pformat

        self.logger.info("Opening config.yaml file ...")
        with open(config, "r") as f:
            parsed_config = yaml.safe_load(f)['era5']
        self.logger.info("Loaded ERA5 configuration:\n%s", pformat(parsed_config))
        return parsed_config

    def get_data(self):
        import os
        import cdsapi
        import logging
        import pandas as pd
        from pprint import pformat

        timerange = pd.Series(pd.date_range(self.tstart, self.tend, freq=self.config['timefreq'], inclusive='left'))
        out_dir = timerange[0].strftime(self.config['out_dir'])
        target = f"{out_dir}/era5_{self.tstart}_{self.tend}.nc"
        if not os.path.exists(os.path.dirname(target)):
            os.makedirs(os.path.dirname(target))
        
        client = cdsapi.Client()
        request = {
            "product_type": [self.config['product_type']],
            "variable": [self.config['variables']],
            "year": [f"{d:02d}" for d in timerange.dt.year.unique()],
            "month": [f"{d:02d}" for d in timerange.dt.month.unique()],
            "day": [f"{d:02d}" for d in timerange.dt.day.unique()],
            "time": [f"{d:02d}:00" for d in timerange.dt.hour.unique()],
            "data_format": self.config['format'],
            "download_format": self.config['download_format'],
            "area": self.config['area']
        }
        self.logger.info("Retrieving ERA5 configuration:\n%s", pformat(request))
        client.retrieve(self.config['dataset'], request).download()

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="""
        ERA5 Downloader. This script requires config.yaml file storing the following parameters:\n
        - dataset-name\n
        - product-type\n
        - variables\n
        - timefreq (max 1hourly)\n
        - area for subsetting\n
        - file format\n
        - download format
        ----------------------------\n
        version 1\n
        created by: @yothunder
        """
    )
    parser.add_argument("tstart", type=str, help="Time start in YYYYMMDD. Example: 20241201", metavar="tstart")
    parser.add_argument("tend", type=str, help="Time end in YYYYMMDD. Example: 20250101", metavar="tend")
    args = parser.parse_args()

    getERA5(args.tstart, args.tend).get_data()