 void sarcomeric_structure(){
            std::vector<vector> myosin_positions;
            box[0] = 16.2;
            box[1] = 12;
            myosin_positions = {{-4.05,-3},{4.05,-3},{-4.05,0},{4.05,0}, {-4.05,3},{4.05,3}};
            for (int i = 0; i < myosin_positions.size(); i++){
                myosin.xs[i][0] = myosin_positions[i][0];
                myosin.xs[i][1] = myosin_positions[i][1];
                myosin.thetas[i] = 0;
            }
            //myosin.n = myosin_positions.size();
            myosin.update_endpoints();
            std::vector<vector> actin_positions;
            actin_positions = {{-6.6,-3.2},{1.5,-3.2},{-6.6,-3},{1.5,-3},{-6.6,-2.8},{1.5,-2.8},{-6.6,-0.2},{1.5,-0.2},{-6.6,0},{1.5,0},{-6.6,0.2},{1.5,0.2},{-6.6,2.8},{1.5,2.8},{-6.6,3},{1.5,3},{-6.6,3.2},{1.5,3.2}};
            for (int i = 0; i < actin_positions.size(); i++){
                actin.xs[i][0] = actin_positions[i][0];
                actin.xs[i][1] = actin_positions[i][1];
                actin.thetas[i] = 0;
            }
            int n = actin_positions.size();
            actin_positions = {{-1.5,-3.2},{6.6,-3.2},{-1.5,-3},{6.6,-3},{-1.5,-2.8},{6.6,-2.8},{-1.5,-0.2},{6.6,-0.2},{-1.5,0},{6.6,0},{-1.5,0.2},{6.6,0.2},{-1.5,2.8},{6.6,2.8},{-1.5,3},{6.6,3},{-1.5,3.2},{6.6,3.2}};
            for (int i = 0; i < actin_positions.size(); i++){
                actin.xs[i+n][0] = actin_positions[i][0];
                actin.xs[i+n][1] = actin_positions[i][1];
                actin.thetas[i+n] = M_PI;
            }
            //actin.n = actin_positions.size()*2;
            actin.update_endpoints();
            update_system();
            //save_state();
        }

        void partial_sarcomeric_structure(){
            std::vector<vector> myosin_positions;
            box[0] = 19.2;
            box[1] = 12;
            myosin_positions = {{-4,-3},{4,-3},{-4,0},{4,0}, {-4,3},{4,3}};
            for (int i = 0; i < myosin_positions.size(); i++){
                myosin.xs[i][0] = myosin_positions[i][0];
                myosin.xs[i][1] = myosin_positions[i][1];
                myosin.thetas[i] = 0;
            }
            //myosin.n = myosin_positions.size();
            myosin.update_endpoints();
            //alpha_actinin.n = alpha_actinin_positions.size();
            std::vector<vector> actin_positions;
            actin_positions = {{1.5,-3.2},{1.5,-3},{1.5,-2.8},{1.5,-0.2},{1.5,0},{1.5,0.2},{1.5,2.8},{1.5,3},{1.5,3.2}};
            for (int i = 0; i < actin_positions.size(); i++){
                actin.xs[i][0] = actin_positions[i][0];
                actin.xs[i][1] = actin_positions[i][1];
                actin.thetas[i] = 0;
            }
            int n = actin_positions.size();
            actin_positions = {{-1.5,-3.2},{-1.5,-3},{-1.5,-2.8},{-1.5,-0.2},{-1.5,0},{-1.5,0.2},{-1.5,2.8},{-1.5,3},{-1.5,3.2}};
            for (int i = 0; i < actin_positions.size(); i++){
                actin.xs[i+n][0] = actin_positions[i][0];
                actin.xs[i+n][1] = actin_positions[i][1];
                actin.thetas[i+n] = M_PI;
            }
            //actin.n = actin_positions.size()*2;
            actin.update_endpoints();
            update_system();
            //save_state();
        }

        void catch_bond(){
            std::vector<vector> myosin_positions;
            box[0] = 19.2;
            box[1] = 12;
            myosin_positions = {{-4,0},{4,0}};
            for (int i = 0; i < myosin_positions.size(); i++){
                myosin.xs[i][0] = myosin_positions[i][0];
                myosin.xs[i][1] = myosin_positions[i][1];
                myosin.thetas[i] = 0;
            }
            //myosin.n = myosin_positions.size();
            myosin.update_endpoints();
            std::vector<vector> alpha_actinin_positions;
            alpha_actinin_positions = {{0,0},{-9.6,0}};
            for (int i = 0; i < alpha_actinin_positions.size(); i++){
                alpha_actinin.xs[i][0] = alpha_actinin_positions[i][0];
                alpha_actinin.xs[i][1] = alpha_actinin_positions[i][1];
            }
            //alpha_actinin.n = alpha_actinin_positions.size();
            std::vector<vector> actin_positions;
            actin_positions = {{1.5,0}};
            for (int i = 0; i < actin_positions.size(); i++){
                actin.xs[i][0] = actin_positions[i][0];
                actin.xs[i][1] = actin_positions[i][1];
                actin.thetas[i] = 0;
            }
            int n = actin_positions.size();
            actin_positions = {{-1.5,0}};
            for (int i = 0; i < actin_positions.size(); i++){
                actin.xs[i+n][0] = actin_positions[i][0];
                actin.xs[i+n][1] = actin_positions[i][1];
                actin.thetas[i+n] = M_PI;
            }
            //actin.n = actin_positions.size()*2;
            actin.update_endpoints();
            update_system();
            //save_state();
        }