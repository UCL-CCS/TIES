import datetime
import pathlib
import subprocess
import zipfile
import logging
import os

from flask import Flask, render_template, request, redirect, send_file

def create_app(storageties, ambertools, load_ties):
    logger = logging.getLogger('waitress')
    logger.setLevel(logging.INFO)

    app = Flask(__name__)
    app.config['SECRET_KEY'] = b'\xe2D\xf9\x90\t\x89\xf9\xaa\x87\x02\xf5$\t\xc4\xd3\xa0'

    # the working directory
    work_dir = pathlib.Path(storageties)

    @app.route('/test')
    def test():
        return 'simple'

    @app.route('/')
    def hello_world():
        return render_template('main.html')

    @app.route("/upload-image", methods=["GET", "POST"])
    def upload_image():
        if request.method == "POST":
            # if request.form['password'] != 'haslo':
            #     return "Wrong Password"

            # net charge
            print('Verifying charge')
            logger.info('Doing the charge -nc testing.')
            not_validated_nc = request.form['net_charge']
            if not_validated_nc.startswith('-'):
                nc_negative = True
                not_validated_nc = not_validated_nc[1:]
            else:
                nc_negative = False
            if not_validated_nc.isdigit():
                net_charge = int(request.form['net_charge'])
            else:
                return 'Net ligand charge is not an integer'
            print('Charge verified')

            # check the files
            if request.files['ligand_ini'].filename == '' or request.files['ligand_fin'].filename == '':
                return 'One of the ligands was not uploaded'

            # create a dedicated directory (date and metadata) for the request where to save the files etc
            session_dir = work_dir / f'{datetime.datetime.now().strftime("%d-%m-%Y-%H:%M:%S:%f")}'
            # todo - save IP address, and other info related in the request (entire request?)
            session_dir.mkdir()

            # save the files
            request.files['ligand_ini'].save(session_dir / request.files['ligand_ini'].filename)
            request.files['ligand_fin'].save(session_dir / request.files['ligand_fin'].filename)

            # redirect?
            # run TIES on the underlying system
            # activate ties env
            loadties = load_ties
            # run it
            source_amber = "export AMBERHOME=/home/ccsadmin/anaconda3/pkgs/ambertools-20.15-py39he1d1048_1 ; "
            command = f'{source_amber} /home/ccsadmin/anaconda3/bin/ties create ' \
                      f'-l {request.files["ligand_ini"].filename} ' \
                      f'{request.files["ligand_fin"].filename} ' \
                      f'-nc {net_charge}'
            print(f'About to try the command: {command}')
            with open(session_dir / 'run.log', 'w') as LOG:
                output = subprocess.run([command], shell=True, stdout=LOG, stderr=LOG, cwd=session_dir)
            # todo apply the /n to be actually not escaped characters?
            print('OUTPUT', output)

            # zip altogether
            os.chdir(session_dir)
            with zipfile.ZipFile(session_dir / 'ties20.zip', 'w') as myzip:
                for i in (session_dir / 'ties20').glob('**/*'):
                    print('rel', i.relative_to((session_dir / 'ties20')))
                    myzip.write(i.relative_to(session_dir))
                # add the main log
                myzip.write(session_dir / 'run.log')
            zipped_output = session_dir / 'ties20.zip'

            return send_file(zipped_output, as_attachment=True,
                             attachment_filename=f'ties20_{request.files["ligand_ini"].filename}_{request.files["ligand_fin"].filename}.zip')

        return render_template("public/upload_image.html")

    return app


if __name__ == '__main__':
    app = create_app('/home/dresio/tiesclients',
                     'source /home/dresio/software/amber18install/amber.sh',
                     'source /home/dresio/software/virtualenvs/tiesdev/bin/activate '
                     )
    app.run()