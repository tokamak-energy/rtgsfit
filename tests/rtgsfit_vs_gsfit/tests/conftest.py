from datetime import datetime
import logging
import os
import subprocess

def pytest_sessionfinish(session, exitstatus):
    if exitstatus == 0:  # All tests passed
        commit_hash = subprocess.check_output(
            ["git", "rev-parse", "HEAD"]
        ).decode().strip()
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

        certificate_text = (
            f"TEST CERTIFICATE\n"
            f"=================\n"
            f"All tests PASSED successfully.\n"
            f"Commit: {commit_hash}\n"
            f"Time: {timestamp}\n"
        )

        with open("test_certificate.txt", "w") as f:
            f.write(certificate_text)

        print("\n" + certificate_text)


def pytest_configure(config):
    # Detect if running in xdist parallel mode
    worker = os.environ.get("PYTEST_XDIST_WORKER", "master")
    logfile = f"pytest_{worker}.log"

    # Basic logging setup
    logging.basicConfig(
        filename=logfile,
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s"
    )

    # Optional: log a startup message
    logging.info(f"Starting tests on worker '{worker}'")
