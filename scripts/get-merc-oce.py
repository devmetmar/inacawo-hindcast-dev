class getMercatorOcean:
    def __init__(self, tstart:str, tend:str):
        self.tstart = tstart
        self.tend = tend
        self.logger = self._setup_logger()
        self.logger.info(f"""
        ============================
        Initialized Mercator Ocean downloader.
        TSTART: {self.tstart}
        TEND  : {self.tend}
        ============================
        """)
        self.config = self.load_yaml()

    def _setup_logger(self):
        import os
        import logging
        from datetime import datetime
        
        logdir  = datetime.strptime(self.tstart, "%Y%m%d").strftime("/home/cawofcst_dev/devmetmar/logs/download/mercatorocean/%Y/%m")
        logname = f"{logdir}/mercatorocean_{self.tstart}_{self.tend}_{datetime.utcnow().strftime('%Y%m%d_%H%M%S')}.log"
        if not os.path.exists(logdir):
            os.makedirs(logdir)
        # logging.basicConfig(
        #     filename=logname,
        #     level=logging.INFO,
        #     format="%(asctime)s - %(levelname)s - %(message)s"
        # )
        logger = logging.getLogger("MERCATOROCEANDownloader")
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
            parsed_config = yaml.safe_load(f)['mercatorocean']
        self.logger.info("Loaded MERCATOR OCEAN configuration:\n%s", pformat(parsed_config))
        return parsed_config

    def get_data(self, **kwargs):
        import os
        import copernicusmarine
        from datetime import datetime
        from pprint import pformat

        if datetime.strptime(self.tstart, "%Y%m%d") < datetime(2021,7,1):
            dataset_id = self.config['dataset_id'][0]
        else:
            dataset_id = self.config['dataset_id'][1]

        out_dir = datetime.strptime(self.tstart, "%Y%m%d").strftime(self.config['out_dir'])
        target = f"mercatorocean_{self.tstart}_{self.tend}.nc"
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

        subset_args = {
            "dataset_id": dataset_id,
            "variables": self.config['variables'],
            "minimum_longitude": self.config['area'][1],
            "maximum_longitude": self.config['area'][3],
            "minimum_latitude": self.config['area'][2],
            "maximum_latitude": self.config['area'][0],
            "start_datetime": self.tstart,
            "end_datetime": self.tend,
            "output_filename": target,
            "output_directory": out_dir,
            "username": self.config['user'],
            "password": self.config['pass']
        }

        subset_args.update(kwargs)
        self.logger.info("Calling copernicusmarine.subset with arguments:\n%s", pformat(subset_args))
        copernicusmarine.subset(**subset_args)

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="""
        MERCATOR OCEAN Downloader. This script requires config.yaml file storing the following parameters:\n
        - dataset-id\n
        - variables\n
        - area for subsetting\n
        - file format\n
        - out_dir\n
        - authentication\n
        ----------------------------\n
        version 1\n
        created by: @yothunder
        """
    )
    parser.add_argument("tstart", type=str, help="Time start in YYYYMMDD. Example: 20241201", metavar="tstart")
    parser.add_argument("tend", type=str, help="Time end in YYYYMMDD. Example: 20250101", metavar="tend")
    args = parser.parse_args()

    getMercatorOcean(args.tstart, args.tend).get_data()