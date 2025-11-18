import logging

logger = logging.getLogger("PLCA")
logger.setLevel(logging.INFO)
ch_format = logging.Formatter('[%(asctime)s] %(levelname)s: %(message)s', "%Y-%m-%d %H:%M:%S")
ch = logging.StreamHandler()
ch.setFormatter(ch_format)
logger.addHandler(ch)

