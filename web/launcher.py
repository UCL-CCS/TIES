import datetime
import pathlib
import subprocess
import zipfile
import logging
import os

from flask import Flask, render_template, request, redirect, send_file

def is_number(n):
    try:
        float(n)
        return True
    except:
        return False

def is_integer(n):
    if is_number(n) and float(n) == int(n):
        return True

    return False

def create_app(storageties, tiesenv):
    logger = logging.getLogger('waitress')
    logger.setLevel(logging.INFO)

    app = Flask(__name__)
    app.config['SECRET_KEY'] = b'\xe2D\xf9\x90\t\x89\xf9\xaa\x87\x02\xf5$\t\xc4\xd3\xa0'

    # the working directory
    work_dir = pathlib.Path(storageties)

    @app.route('/test')
    def test():
        return 'simple test'

    @app.route("/", methods=["GET", "POST"])
    def upload_image():
        if request.method == "POST":
            if request.form['password'] != 'xQUqUzAO':
                return "Wrong Password"

            # net charge
            print('Verifying charge')
            logger.info('Doing the charge -nc testing.')
            if is_integer(request.form['net_charge']):
                net_charge = int(request.form['net_charge'])
            else:
                return 'Net ligand charge is not an integer'
            print(f'Net Ligand Charge verified: {net_charge}')

            # q tol
            if is_number(request.form['q_tol']):
                q_tol = float(request.form['q_tol'])
            else:
                return 'Q tol is not a number'

            # Net q tol
            if is_number(request.form['net_q_tol']):
                net_q_tol = float(request.form['net_q_tol'])
            else:
                return 'Net q tol is not a number'

            # check the files
            if request.files['ligand_ini'].filename == '' or request.files['ligand_fin'].filename == '':
                return 'One of the ligands was not uploaded'

            # create a dedicated directory (date and metadata) for the request where to save the files etc
            session_dir = work_dir / f'{datetime.datetime.now().strftime("%d-%m-%Y-%H:%M:%S:%f")}'
            session_dir.mkdir()
            # todo - save IP address, and other info related in the request (entire request?)

            # save the files
            request.files['ligand_ini'].save(session_dir / request.files['ligand_ini'].filename)
            request.files['ligand_fin'].save(session_dir / request.files['ligand_fin'].filename)

            command = f'{tiesenv} ties create ' \
                      f'--ligands {request.files["ligand_ini"].filename} {request.files["ligand_fin"].filename} ' \
                      f'--ligand-net-charge {net_charge} ' \
                      f'--q-pair-tolerance {q_tol} ' \
                      f'--q-net-tolerance {net_q_tol} '
            print(f'About to run the command: {command}')
            with open(session_dir / 'ties20.log', 'w') as LOG:
                subprocess.run([command], shell=True, stdout=LOG, stderr=LOG, cwd=session_dir)

            # zip the output
            os.chdir(session_dir)
            with zipfile.ZipFile(session_dir / 'ties20.zip', 'w') as myzip:
                for i in (session_dir / 'ties20').glob('**/*'):
                    print('rel', i.relative_to((session_dir / 'ties20')))
                    myzip.write(i.relative_to(session_dir))
                # add the main log
                myzip.write('ties20.log')
            zipped_output = session_dir / 'ties20.zip'

            # send the file back
            return send_file(zipped_output, as_attachment=True,
                             attachment_filename=f'ties20_{request.files["ligand_ini"].filename}_{request.files["ligand_fin"].filename}.zip')

        return render_template("main.html")

    return app


if __name__ == '__main__':
    app = create_app('/home/dresio/tiesclients',
                     'source /home/dresio/software/amber18install/amber.sh ; ' +
                     'source /home/dresio/software/virtualenvs/tiesdev/bin/activate ; '
                     )
    app.run()