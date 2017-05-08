REM celery -A cueBeamCore2 -b "amqp://guest:guest@192.168.0.43:5672/" worker --loglevel info
celery -A cueBeamCore2 -b "amqp://guest:guest@192.168.0.43:5672/" worker
