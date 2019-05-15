#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from PRIMeval_remote import primer_run, upload_results, update_status
from multiprocessing import Process
from rq import Queue, Worker, Connection, get_failed_queue
from rq.job import Job
from rq.registry import FinishedJobRegistry
import shutil
import random
import urllib3
import tarfile
import redis
import time
import base64
import os

__author__ = "Rick Conzemius"
__copyright__ = "Copyright 2019, Rick Conzemius"
__license__ = "MIT"
__maintainer__ = "Rick Conzemius"
__version__ = "1.0"
__email__ = "rick.conzemius.fl@ait.ac.at"
__status__ = "Production"

key = "publickey"
download_url = "http://download_url/" + key
download_folder = os.getcwd() + "/downloads/"
run_folder = os.getcwd() + "/runs/"

redisClient = redis.Redis(host='localhost', port=6379, db=0)
RunQueue = Queue('RunQueue',connection=redisClient)
UploadQueue = Queue('UploadQueue',connection=redisClient)

def download_runs():
    tmp_file = "download_" + str(random.randint(100000, 999999))
    tmp_path = download_folder + tmp_file
    http = urllib3.PoolManager()
    try:
        with http.request('GET', download_url, preload_content=False) as r, open(tmp_path, 'wb') as out_file:
            shutil.copyfileobj(r, out_file)
    except:
        print("Download runs connection error.")
        return

    if r.status == 200:
        file_name = r.headers['content-disposition'].split("=", 1)[1]
        run_id = file_name.split(".", 1)[0]
        current_run_folder = run_folder + str(run_id) + "/"
        update_status(run_id, "Download", "Download", 0, 1)
        update_status(run_id, "Download", "Extraction", 1, 0)

        try:
            archive = tarfile.open(tmp_path)
            archive.extractall(run_folder)
            archive.close()
            tar_error = False
        except tarfile.ReadError:
            tar_error = True


        if tar_error == False:
            for seqtype in ['primers', 'probes', 'contigs']:
                if not os.path.exists(current_run_folder + "input/" + seqtype):
                    os.makedirs(current_run_folder + "input/" + seqtype)
            update_status(run_id, "Download", "Extraction", 0, 1)

            RunQueue.enqueue(primer_run, args = (run_id,), job_id = str(run_id), timeout=7200)
            update_status(run_id, "PRIMEval", "PreparingFolders", 1, 0)
        else:
            update_status(run_id, "Download", "Extraction", 0, 0)
            if os.path.exists(current_run_folder):
                os.remove(current_run_folder)

    if os.path.exists(tmp_path):
        os.remove(tmp_path)

def failed_jobs_by_queue(queue):
    failed = []
    fq = get_failed_queue(redisClient)
    for id in fq.get_job_ids():
        try:
           if Job.fetch(str(id), connection = redisClient).origin == queue.name:
                failed.append(str(id))
        except:
            continue
    return failed

def check_run_queue():
    # Failed jobs
    for id in failed_jobs_by_queue(RunQueue):
        update_status(id, "PRIMEval", "Run", 0, 0)
        if os.path.exists(run_folder + str(id)):
            shutil.rmtree(run_folder + str(id))
        try:
            error_message = Job.fetch(str(id), connection=redisClient).exc_info
            error_message = error_message[-1000:]
            error_message = base64.urlsafe_b64encode(error_message.encode()).decode()
            Job.fetch(str(id), connection=redisClient).delete()
            update_status(id, "PRIMEval", "Run", 0, 0, error_message)
        except:
            return

    # Finished jobs
    for id in FinishedJobRegistry('RunQueue', connection=redisClient).get_job_ids():
        update_status(id, "PRIMEval", "Run", 0, 1)
        try:
            Job.fetch(str(id), connection=redisClient).delete()
        except:
            return
        UploadQueue.enqueue(upload_results, args = (id,), job_id = str(id), timeout = 1200)
        update_status(id, "Upload", "Queued", 1, 0)

def check_upload_queue():
    # Failed jobs
    for id in failed_jobs_by_queue(UploadQueue):
        update_status(id, "Upload", "Upload", 0, 0)
        if os.path.exists(run_folder + str(id)):
            shutil.rmtree(run_folder + str(id))
        try:
            Job.fetch(str(id), connection=redisClient).delete()
        except:
            return
    # Finished jobs
    for id in FinishedJobRegistry('UploadQueue', connection=redisClient).get_job_ids():
        if os.path.exists(run_folder + str(id)):
            shutil.rmtree(run_folder + str(id))
        try:
            Job.fetch(str(id), connection=redisClient).delete()
        except:
            return

def manage_queues():
    while True:
        print("Checking queues.")
        check_run_queue()
        check_upload_queue()
        time.sleep(5)

def manage_downloads():
    while True:
        print("Downloading run.")
        download_runs()
        time.sleep(10)

def runqueue_worker():
    with Connection(redisClient):
        worker = Worker(map(Queue, ["RunQueue"]))
        worker.work()

def uploadqueue_worker():
    with Connection(redisClient):
        worker = Worker(map(Queue, ["UploadQueue"]))
        worker.work()

try:
    p1 = Process(target=runqueue_worker)
    p2 = Process(target=uploadqueue_worker)
    p3 = Process(target=manage_queues)
    p4 = Process(target=manage_downloads)
    p1.start()
    p2.start()
    time.sleep(2)
    p3.start()
    p4.start()
except KeyboardInterrupt:
    p1.terminate()
    p2.terminate()
    p3.terminate()
    p4.terminate()
