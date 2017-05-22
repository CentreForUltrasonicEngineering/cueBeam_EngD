REM celery -A cueBeamCore3 -b "amqp://guest:guest@EEE-mimir.ds.strath.ac.uk:5672/" worker --loglevel info
celery -A cueBeamCore3 -b "amqp://guest:guest@EEE-mimir.ds.strath.ac.uk:5672/" worker
